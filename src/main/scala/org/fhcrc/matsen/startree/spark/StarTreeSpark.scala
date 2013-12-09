package org.fhcrc.matsen.startree.spark

// Start
import java.io._
import scala.collection.JavaConverters._

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMRecord

import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.evolution.alignment.Alignment;
import org.fhcrc.matsen.startree._;

import com.google.common.base.Preconditions;

object StarTreeSpark {
  def main(args: Array[String]) {
    // spark path, json fasta, sam
    if(args.length != 4) {
      System.err.format("USAGE: StarTreeSpark master_path json fasta bam\n");
      System.exit(1);
    }

    val masterPath = args(0)
    val jsonPath = args(1)
    val fastaPath = args(2)
    val bamPath = args(3)

    val reader = new SAMFileReader(new File(bamPath))
    System.err.format("Loading JSON from %s\n", jsonPath);
    val jsonReader = new BufferedReader(new FileReader(jsonPath));
    val modelRates = HKYModelParser.substitutionModel(jsonReader).asScala;
    jsonReader.close()

    System.setProperty("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
    System.setProperty("spark.kryo.registrator", "org.fhcrc.matsen.startree.spark.StarTreeKryoRegistrator");
    System.setProperty("spark.kryoserializer.buffer.mb", "64");

    val sc = new SparkContext(masterPath, "StarTreeRenaissance")
    //val sc = new SparkContext("local", "Simple App", "YOUR_SPARK_HOME",
    //  List("target/scala-2.9.3/simple-project_2.9.3-1.0.jar"))

    val references = SAMUtils.readAllFasta(new File(fastaPath))

    val alignments = reader.asScala.map(r => {
      val ref = references.get(r.getReferenceName());
      Preconditions.checkNotNull(ref, "No reference for %s", r.getReferenceName());
      SAMBEASTUtils.alignmentOfRecord(r, ref)
    }).toList

    val byReference = sc.parallelize(alignments)
      .keyBy(a => a.getTaxon(0).getId)
      .map { case (ref, a) => {
        val model = modelRates.map(hr => hr.getModel).asJava
        val rates = modelRates.map(hr => hr.getSiteRateModel).asJava
        (ref, StarTreeRenaissance.calculate(a, model, rates))
      } }
      .reduceByKey(_.plus(_)).collect
  }
}
