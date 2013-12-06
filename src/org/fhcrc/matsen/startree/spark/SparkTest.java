package org.fhcrc.matsen.startree.spark;

import beagle.BeagleInfo;
import com.google.common.base.Function;
import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.alignment.Alignment;
import dr.inference.model.Parameter;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function2;
import org.fhcrc.matsen.startree.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by cmccoy on 12/6/13.
 */
public class SparkTest {
    public static void main(String... args) throws Exception {
        // spark path, json fasta, sam
        if(args.length != 4) {
            System.err.format("USAGE: SparkTest master_path json fasta bam\n");
            System.exit(1);
        }

        final String masterPath = args[0],
                jsonPath = args[1],
                fastaPath = args[2],
                bamPath = args[3];

        SAMFileReader reader = new SAMFileReader(new File(bamPath));

        BeagleInfo.getVersion();

        System.err.format("Loading JSON from %s\n", jsonPath);
        BufferedReader jsonReader = new BufferedReader(new FileReader(jsonPath));
        final List<HKYModelParser.HKYAndRate> mRates = HKYModelParser.substitutionModel(jsonReader);

        System.setProperty("spark.serializer", "org.apache.spark.serializer.KryoSerializer");
        JavaSparkContext ctx = new JavaSparkContext(masterPath, "StarTreeRenaissance");
        final Map<String, byte[]> references = SAMUtils.readAllFasta(new File(fastaPath));

        List<Alignment> alignments = Lists.newArrayList(Iterables.transform(reader, new Function<SAMRecord, Alignment>() {
            @Override
            public Alignment apply(SAMRecord samRecord) {
                final byte[] ref = references.get(samRecord.getReferenceName());
                Preconditions.checkNotNull(ref, "no reference for %s", samRecord.getReferenceName());
                return SAMBEASTUtils.alignmentOfRecord(samRecord, ref);
            }
        }));

        final TwoTaxonResult result = ctx.parallelize(alignments).map(new org.apache.spark.api.java.function.Function<Alignment, TwoTaxonResult>() {
            @Override
            public TwoTaxonResult call(final Alignment a) throws Exception {
                final List<HKY> models = new ArrayList<HKY>();
                final List<SiteRateModel> rates = new ArrayList<SiteRateModel>();
                for (HKYModelParser.HKYAndRate hr : mRates) {
                    models.add(hr.getModel());
                    rates.add(hr.getSiteRateModel());
                }

                System.err.println(System.getProperty("java.library.path"));

                BeagleInfo.printResourceList();

                return StarTreeRenaissance.calculate(a, models, rates);
            }
        }).reduce(new Function2<TwoTaxonResult, TwoTaxonResult, TwoTaxonResult>() {
            @Override
            public TwoTaxonResult call(TwoTaxonResult a, TwoTaxonResult b) throws Exception {
                return a.plus(b);
            }
        });
    }
}
