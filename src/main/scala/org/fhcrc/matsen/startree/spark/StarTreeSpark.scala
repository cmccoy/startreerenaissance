package org.fhcrc.matsen.startree.spark

// Start

import java.io._
import java.util.logging._
import scala.collection.JavaConverters._

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import net.sf.samtools.SAMFileReader

import org.fhcrc.matsen.startree._

import com.google.common.base.Preconditions
import com.amazonaws.services.s3.AmazonS3Client

case class StarTreeSparkConfig(parallelism: Int = 12,
                               prefix: String = "",
                               smooth: Boolean = true,
                               sample: Boolean = false,
                               masterPath: String = "",
                               sparkHome: String = "/root/spark",
                               jarPath: String = "/root/startreerenaissance/target/scala-2.9.3/startreerenaissance-assembly-0.1.jar",
                               consolidateFiles: Boolean = true,
                               executorMemory: String = "8g",
                               bucket: Option[String] = None,
                               jsonPath: File = new File("."),
                               fastaPath: File = new File("."),
                               bamPath: File = new File("."))

object StarTreeSpark {
  val appName = "StarTreeSpark"
  val logger = Logger.getLogger("org.fhcrc.matsen.startree.spark")

  val parser = new scopt.OptionParser[StarTreeSparkConfig](appName) {
    head(appName, "0.1")
    opt[Int]('j', "parallelism") action {
      (x, c) => c.copy(parallelism = x)
    } text ("Parallelism level")
    opt[String]('p', "prefix") action {
      (x, c) => c.copy(prefix = x)
    } text ("Prefix to add to each output file")
    opt[Unit]('n', "no-smooth") action {
      (_, c) => c.copy(smooth = false)
    } text ("Do *not* smooth parameter estimates using empirical bayes")
    opt[Unit]("no-consolidate-files") action {
      (_, c) => c.copy(consolidateFiles = false)
    } text ("Do *not* consolidate files")
    opt[Unit]('n', "sample") action {
      (_, c) => c.copy(sample = true)
    } text ("Sample rates from poisson-gamma, rather than just using the mean")
    opt[String]('b', "bucket") action {
      (x, c) => c.copy(bucket = Some(x))
    } text ("Upload to bucket")
    opt[String]('s', "spark-home") action {
      (x, c) => c.copy(sparkHome = x)
    } text ("Spark home")
    opt[String]('j', "jar-path") action {
      (x, c) => c.copy(jarPath = x)
    } text ("Path to the application JAR")
    opt[String]("executor-memory") action {
      (x, c) => c.copy(executorMemory = x)
    } text ("Executor memory")
    arg[String]("<master_path>") required() action {
      (x, c) => c.copy(masterPath = x)
    } text ("path to SPARK master")
    arg[File]("<json>") required() action {
      (x, c) => c.copy(jsonPath = x)
    } text ("path to JSON model specification")
    arg[File]("<fasta>") required() action {
      (x, c) => c.copy(fastaPath = x)
    } text ("path to reference FASTA file")
    arg[File]("<bam>") required() action {
      (x, c) => c.copy(bamPath = x)
    } text ("path to BAM file with aligned reads.")
    help("help") text ("prints this usage text")
  }

  def run(config: StarTreeSparkConfig) {
    val reader = new SAMFileReader(config.bamPath)
    System.err.format("Loading JSON from %s\n", config.jsonPath);
    val jsonReader = new BufferedReader(new FileReader(config.jsonPath));
    val modelRates = HKYModelParser.substitutionModel(jsonReader).asScala;
    jsonReader.close()

    System.setProperty("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
    System.setProperty("spark.kryo.registrator", "org.fhcrc.matsen.startree.spark.StarTreeKryoRegistrator")
    System.setProperty("spark.kryoserializer.buffer.mb", "256")
    System.setProperty("spark.executor.memory", config.executorMemory)
    System.setProperty("spark.akka.frameSize", "512")
    if (config.consolidateFiles)
      System.setProperty("spark.shuffle.consolidateFiles", "true")

    val sc = config.masterPath match {
      case mp if mp.startsWith("local") =>
        new SparkContext(mp, appName)
      case mp =>
        new SparkContext(mp, appName,
          config.sparkHome,
          Seq(config.jarPath))
    }

    val references = SAMUtils.readAllFasta(config.fastaPath)

    val alignments = sc.parallelize(reader.asScala.map(
      r => {
        val ref = references.get(r.getReferenceName());
        Preconditions.checkNotNull(ref, "No reference for %s", r.getReferenceName());
        SAMBEASTUtils.alignmentOfRecord(r, ref)
      }).toList, config.parallelism).keyBy(_.getTaxon(0).getId)

    val result = alignments.mapValues(a => {
      val model = modelRates.map(hr => hr.getModel).asJava
      val rates = modelRates.map(hr => hr.getSiteRateModel).asJava
      StarTreeRenaissance.calculate(a, model, rates)
    }).reduceByKey(_.plus(_), config.parallelism)

    val bcastConfig = sc.broadcast(config)

    val stResults = result.map {
      case (refName, v) =>
        (refName,
         new StarTreeRenaissanceResult(refName, v, bcastConfig.value.sample, StarTreeRenaissance.CHAIN_LENGTH / 10))
    }.cache

    // Save each target to S3 if a bucket is given
    if (!config.bucket.isEmpty) {
      println("Uploading to S3 bucket " + config.bucket.get)
      val cred = sc.broadcast(new com.amazonaws.auth.DefaultAWSCredentialsProviderChain().getCredentials)

      stResults foreach {
        case (refName, result) => {
          val jsonBase = config.prefix + refName.replaceAll("[*/]+", "_")
          val jsonName = jsonBase + ".json.gz"
          val tmpFile = File.createTempFile(jsonBase, ".json.gz")

          val jsonStream = new java.util.zip.GZIPOutputStream(new FileOutputStream(tmpFile))
          val jsonWriter = new PrintStream(jsonStream)
          val gson = org.fhcrc.matsen.startree.gson.getGsonBuilder.create
          gson.toJson(result, jsonWriter)
          jsonWriter.close()

          // errors on missing credentials
          val s3Client = new AmazonS3Client(cred.value)

          try {
            val r = s3Client.putObject(bcastConfig.value.bucket.get, jsonName, tmpFile)
            println(r)
          } finally {
            val wasDeleted = tmpFile.delete()
            require(wasDeleted, "failed to delete temp file")
          }
        }
      }
    }

    val resultLocal = stResults.collect

    // Stop the Spark context - remaining work is local.
    sc.stop()

    resultLocal foreach {
      case (refName, result) => {
        val sanitized = refName.replaceAll("[*/]+", "_")
        val outName = config.prefix + sanitized + ".log"
        val writer = new PrintStream(new File(outName))
        if (config.smooth)
          result.getSmoothed.print(writer, true)
        else
          result.getUnsmoothed.print(writer, true)
        writer.close()

        val jsonName = config.prefix + sanitized + ".json.gz"
        val jsonStream = new java.util.zip.GZIPOutputStream(new java.io.FileOutputStream(jsonName))
        val jsonWriter = new PrintStream(jsonStream)
        val gson = org.fhcrc.matsen.startree.gson.getGsonBuilder.create
        gson.toJson(result, jsonWriter)
        jsonWriter.close()
      }
    }
  }

  def main(args: Array[String]) {
    // spark path, json fasta, sam
    parser.parse(args, StarTreeSparkConfig()) map {
      config =>
      // do stuff
        run(config)
    } getOrElse {
      // arguments are bad, error message will have been displayed
      System.exit(1)
    }
  }

}
