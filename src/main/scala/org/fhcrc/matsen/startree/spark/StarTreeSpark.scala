//package org.fhcrc.matsen.startree.spark;
//
//import beagle.BeagleInfo;
//import com.google.common.base.Function;
//import com.google.common.base.Preconditions;
//import com.google.common.collect.Iterables;
//import com.google.common.collect.Lists;
//import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
//import dr.app.beagle.evomodel.substmodel.HKY;
//import dr.evolution.alignment.Alignment;
//import net.sf.samtools.SAMFileReader;
//import net.sf.samtools.SAMRecord;
//import org.apache.spark.api.java.JavaPairRDD;
//import org.apache.spark.api.java.JavaRDD;
//import org.apache.spark.api.java.JavaSparkContext;
//import org.apache.spark.api.java.function.DoubleFunction;
//import org.apache.spark.api.java.function.Function2;
//import org.fhcrc.matsen.startree.*;
//import scala.Tuple2;
//
//import java.io.BufferedReader;
//import java.io.File;
//import java.io.FileReader;
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Map;
//
///**
// * Created by cmccoy on 12/6/13.
// */
//public class StarTreeSpark {
//    public static void main(String... args) throws Exception {
//        // spark path, json fasta, sam
//        if(args.length != 4) {
//            System.err.format("USAGE: StarTreeSpark master_path json fasta bam\n");
//            System.exit(1);
//        }
//
//        final String masterPath = args[0],
//                jsonPath = args[1],
//                fastaPath = args[2],
//                bamPath = args[3];
//
//        SAMFileReader reader = new SAMFileReader(new File(bamPath));
//
//        //BeagleInfo.getVersion();
//
//        System.err.format("Loading JSON from %s\n", jsonPath);
//        BufferedReader jsonReader = new BufferedReader(new FileReader(jsonPath));
//        final List<HKYModelParser.HKYAndRate> modelRates = HKYModelParser.substitutionModel(jsonReader);
//
//        System.setProperty("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
//        System.setProperty("spark.kryo.registrator", "org.fhcrc.matsen.startree.spark.StarTreeKryoRegistrator");
//        System.setProperty("spark.kryoserializer.buffer.mb", "64");
//
//        JavaSparkContext ctx = new JavaSparkContext(masterPath, "StarTreeRenaissance");
//        final Map<String, byte[]> references = SAMUtils.readAllFasta(new File(fastaPath));
//
//        JavaRDD<Alignment> alignments = ctx.parallelize(Lists.newArrayList(Iterables.transform(reader, new Function<SAMRecord, Alignment>() {
//            @Override
//            public Alignment apply(SAMRecord samRecord) {
//                final byte[] ref = references.get(samRecord.getReferenceName());
//                Preconditions.checkNotNull(ref, "no reference for %s", samRecord.getReferenceName());
//                return SAMBEASTUtils.alignmentOfRecord(samRecord, ref);
//            }
//        })));
//
//        //JavaPairRDD<String, Alignment> alignmentByRef = JavaPairRDD.fromRDD(alignments.keyBy(new org.apache.spark.api.java.function.Function<Alignment, String>() {
//            //@Override
//            //public String call(Alignment a) throws Exception {
//                //return a.getTaxon(0).getId();
//            //}
//        //}));
//
//        //ctx.parallelize(alignments).keyBy(new org.apache.spark.api.java.function.Function<Alignment, String>() {
//            //@Override
//            //public String call(Alignment a) throws Exception {
//                //return a.getTaxon(0).getId();
//            //}
//        //}).map(new org.apache.spark.api.java.function.Function<Tuple2<String, Alignment>, Tuple2<String, TwoTaxonResult>>() {
//            //@Override
//            //public Tuple2<String, TwoTaxonResult> call(Tuple2<String, Alignment> tup) throws Exception {
//                //final List<HKY> models = new ArrayList<HKY>();
//                //final List<SiteRateModel> rates = new ArrayList<SiteRateModel>();
//                //for (HKYModelParser.HKYAndRate hr : modelRates) {
//                    //models.add(hr.getModel());
//                    //rates.add(hr.getSiteRateModel());
//                //}
//
//                //return new Tuple2<String, TwoTaxonResult>(tup._1(), StarTreeRenaissance.calculate(tup._2(), models, rates));
//            //}
//        //});
//        final TwoTaxonResult result = alignments.map(new org.apache.spark.api.java.function.Function<Alignment, TwoTaxonResult>() {
//            @Override
//            public TwoTaxonResult call(final Alignment a) throws Exception {
//                final List<HKY> models = new ArrayList<HKY>();
//                final List<SiteRateModel> rates = new ArrayList<SiteRateModel>();
//                for (HKYModelParser.HKYAndRate hr : modelRates) {
//                    models.add(hr.getModel());
//                    rates.add(hr.getSiteRateModel());
//                }
//
//                return StarTreeRenaissance.calculate(a, models, rates);
//            }
//        }).reduce(new Function2<TwoTaxonResult, TwoTaxonResult, TwoTaxonResult>() {
//            @Override
//            public TwoTaxonResult call(TwoTaxonResult a, TwoTaxonResult b) throws Exception {
//                return a.plus(b);
//            }
//        });
//
//        //result.print(System.out);
//    }
//}

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
    val modelRates = HKYModelParser.substitutionModel(jsonReader);
    jsonReader.close()

    System.setProperty("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
    System.setProperty("spark.kryo.registrator", "org.fhcrc.matsen.startree.spark.StarTreeKryoRegistrator");
    System.setProperty("spark.kryoserializer.buffer.mb", "64");

    val sc = new SparkContext(masterPath, "StarTreeRenaissance")
    //val sc = new SparkContext("local", "Simple App", "YOUR_SPARK_HOME",
    //  List("target/scala-2.9.3/simple-project_2.9.3-1.0.jar"))

    val references = SAMUtils.readAllFasta(new File(fastaPath))
    
    val alignments = sc.parallelize(reader.map(r => {
      val ref = references.get(r.getReferenceName());
      Preconditions.checkNotNull(ref, "No reference for %s", r.getReferenceName());
      SAMBEASTUtils.alignmentOfRecord(r, ref)
    }))
  }
}
