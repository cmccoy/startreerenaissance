package org.fhcrc.matsen.startree.spark

// Start
import java.io._
import java.util.logging._
import scala.collection.JavaConverters._

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMRecord

import dr.app.beagle.evomodel.sitemodel.SiteRateModel
import dr.app.beagle.evomodel.substmodel.HKY
import dr.evolution.alignment.Alignment
import org.fhcrc.matsen.startree._
import org.fhcrc.matsen.startree.gson._

import com.google.common.base.Preconditions
import com.google.gson.GsonBuilder

case class Config(parallelism: Int = 12,
                  prefix: String = "",
                  smooth: Boolean = true,
                  sample: Boolean = false,
                  masterPath: String = "",
                  executorMemory: String = "16g",
                  jsonPath: File = new File("."),
                  fastaPath: File = new File("."),
                  bamPath: File = new File("."))

object StarTreeSpark {
  val appName = "StarTreeSpark"
  val logger = Logger.getLogger("org.fhcrc.matsen.startree.spark")

  val parser = new scopt.OptionParser[Config](appName) {
    head(appName, "0.1")
    opt[Int]('j', "parallelism") action { (x, c) => c.copy(parallelism = x) } text("Parallelism level")
    opt[String]('p', "prefix") action { (x, c) => c.copy(prefix = x) } text("Prefix to add to each output file")
    opt[Unit]('n', "no-smooth") action { (_, c) => c.copy(smooth = false) } text("Do *not* smooth parameter estimates using empirical bayes")
    opt[Unit]('n', "sample") action { (_, c) => c.copy(sample = true) } text("Sample rates from poisson-gamma, rather than just using the mean")
    opt[String]("executor-memory") action { (x, c) => c.copy(executorMemory = x) } text("Prefix to add to each output file")
    arg[String]("<master_path>") required() action { (x, c) => c.copy(masterPath = x) } text("path to SPARK master")
    arg[File]("<json>") required() action { (x, c) => c.copy(jsonPath = x) } text("path to JSON model specification")
    arg[File]("<fasta>") required() action { (x, c) => c.copy(fastaPath = x) } text("path to reference FASTA file")
    arg[File]("<bam>") required() action { (x, c) => c.copy(bamPath = x) } text("path to BAM file with aligned reads.")
    help("help") text("prints this usage text")
  }

  def run(config: Config) {
    val reader = new SAMFileReader(config.bamPath)
    System.err.format("Loading JSON from %s\n", config.jsonPath);
    val jsonReader = new BufferedReader(new FileReader(config.jsonPath));
    val modelRates = HKYModelParser.substitutionModel(jsonReader).asScala;
    jsonReader.close()

    System.setProperty("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
    System.setProperty("spark.kryo.registrator", "org.fhcrc.matsen.startree.spark.StarTreeKryoRegistrator");
    System.setProperty("spark.kryoserializer.buffer.mb", "256");
    System.setProperty("spark.executor.memory", config.executorMemory);
    System.setProperty("spark.akka.frameSize", "512");

    val sc = config.masterPath match {
      case mp if mp.startsWith("local") =>
        new SparkContext(mp, appName)
      case mp =>
        new SparkContext(mp, appName,
                         System.getenv("SPARK_HOME"),
                         Seq(System.getenv("STARTREE_JAR")))
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
      }).reduceByKey(_.plus(_), config.parallelism).collect

    // Stop the Spark context - remaining work is local.
    sc.stop()

    result.foreach {
        case (refName, v) => {
          val outName = config.prefix + refName.replaceAll("\\*", "_") + ".log"
          val writer = new PrintStream(new File(outName))
          lazy val smoothed = v.getSmoothed(config.sample)
          if(config.smooth)
            smoothed.print(writer, true)
          else
            v.print(writer, true)
          writer.close()

          val jsonName = config.prefix + refName.replaceAll("\\*", "_") + ".json.gz"
          val jsonStream = new java.util.zip.GZIPOutputStream(new java.io.FileOutputStream(jsonName))
          val jsonWriter = new PrintStream(jsonStream)
          val gson = new GsonBuilder()
            .registerTypeAdapter(BlockRealMatrixSerializer.serializedType, new BlockRealMatrixSerializer)
            .serializeSpecialFloatingPointValues()
            .create
          val result = Map("unsmoothed" -> v, "smoothed" -> smoothed).asJava
          gson.toJson(result, jsonWriter)
          jsonWriter.close
        }
      }
  }

  def main(args: Array[String]) {
    // spark path, json fasta, sam
    parser.parse(args, Config()) map { config =>
      // do stuff
      run(config)
    } getOrElse {
      // arguments are bad, error message will have been displayed
      System.exit(1)
    }
  }

}
