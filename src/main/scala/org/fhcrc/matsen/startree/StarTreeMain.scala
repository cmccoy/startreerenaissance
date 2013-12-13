package org.fhcrc.matsen.startree

import java.io._
import java.util.logging._
import scala.collection.JavaConverters._

import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMRecord

import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.evolution.alignment.Alignment;

import com.google.common.base.Preconditions;

case class Config(smooth: Boolean = true,
                  jsonPath: File = new File("."),
                  fastaPath: File = new File("."),
                  bamPath: File = new File("."),
                  outputPath: File = new File("."))

object StarTreeMain {
  val appName = "StarTreeMain"
  val logger = Logger.getLogger("org.fhcrc.matsen.startree")

  val parser = new scopt.OptionParser[Config](appName) {
    head(appName, "0.1")
    opt[Unit]('n', "no-smooth") action { (_, c) => c.copy(smooth = false) } text("Do *not* smooth parameter estimates using empirical bayes")
    arg[File]("<json>") required() action { (x, c) => c.copy(jsonPath = x) } text("path to JSON model specification")
    arg[File]("<fasta>") required() action { (x, c) => c.copy(fastaPath = x) } text("path to reference FASTA file")
    arg[File]("<bam>") required() action { (x, c) => c.copy(bamPath = x) } text("path to BAM file with aligned reads.")
    arg[File]("<output>") required() action { (x, c) => c.copy(outputPath = x) } text("path to output file with aligned reads.")
    help("help") text("prints this usage text")
  }

  def run(config: Config) {
    val reader = new SAMFileReader(config.bamPath)
    System.err.format("Loading JSON from %s\n", config.jsonPath);
    val jsonReader = new BufferedReader(new FileReader(config.jsonPath));
    val modelRates = HKYModelParser.substitutionModel(jsonReader).asScala;
    jsonReader.close()

    val references = SAMUtils.readAllFasta(config.fastaPath)

    val alignments = reader.asScala.map(
      r => {
        val ref = references.get(r.getReferenceName());
        Preconditions.checkNotNull(ref, "No reference for %s", r.getReferenceName());
        SAMBEASTUtils.alignmentOfRecord(r, ref)
      }).toList

    val v = alignments.map(a => {
        val model = modelRates.map(hr => hr.getModel).asJava
        val rates = modelRates.map(hr => hr.getSiteRateModel).asJava
        System.err.println("Working on " + a.getTaxon(1).getId())
        StarTreeRenaissance.calculate(a, model, rates)
    }).reduce(_.plus(_))

    val writer = new PrintStream(config.outputPath)
    if(config.smooth)
      v.getSmoothed.print(writer, true)
    else
      v.print(writer, true)
    writer.close()
  }

  def main(args: Array[String]) {
    parser.parse(args, Config()) map { config =>
      // do stuff
      run(config)
    } getOrElse {
      // arguments are bad, error message will have been displayed
      System.exit(1)
    }
  }

}
