package org.fhcrc.matsen.startree;

import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by cmccoy on 12/5/13.
 */
public class StarTreeRenaissanceTest {
    @Before
    public void setUp() throws Exception {

    }

    @After
    public void tearDown() throws Exception {

    }

    //@org.junit.Ignore("Just runs a sampler")
    @Test
    public void testCalculate() throws Exception {
        List<HKY> hkys = new ArrayList<HKY>();
        List<SiteRateModel> rates = new ArrayList<SiteRateModel>();

        for(int i = 0; i < 3; i++) {
            hkys.add(new HKY(1.0, new FrequencyModel(Nucleotides.INSTANCE, new double[]{0.25, 0.25, 0.25, 0.25})));
            final GammaSiteRateModel r = new GammaSiteRateModel(String.format("rate%d", i));
            rates.add(r);
        }

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.addSequence(new Sequence(new Taxon("ref"), "GTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGA"));
        alignment.addSequence(new Sequence(new Taxon("qry"), "GTCACCATGACCACAGACACATCCACGAGCACAGCCCACCTGGAACTGAAGAGCCTGAGATCTGACGACACGGCCGTGTATTTCTGTGCGCGA"));

        TwoTaxonResult r = StarTreeRenaissance.calculate(alignment, hkys, rates, 20000, 200);
        //System.err.format("uN:\n%s\nuS:\n%s\ncN:\n%s\ncS:\n%s\n",
                //r.getUnconditionalNonsynonymous(),
                //r.getUnconditionalSynonymous(),
                //r.getConditionalNonsynonymous(),
                //r.getConditionalSynonymous());
    }
}
