package org.fhcrc.matsen.startree;

import dr.app.beagle.evomodel.substmodel.HKY;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.List;

/**
 * Created by cmccoy on 12/5/13.
 */
public class HKYModelParserTestCase {
    Reader reader = null;
    @Before
    public void setUp() throws Exception {
        reader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("sample_model2.json")));
    }
    @After
    public void tearDown() throws Exception {
        if(reader != null)
            reader.close();
    }

    @Test
    public void testSubstitutionModel() throws Exception {
        List<HKYModelParser.HKYAndRate> result = HKYModelParser.substitutionModel(this.reader);

        assertEquals(3, result.size());
        assertEquals(4, result.get(0).getModel().getFrequencyModel().getFrequencyCount());
        assertEquals(1.4436320149021, result.get(1).getModel().getKappa(), 1e-6);
        assertEquals(0.3199744637909102, result.get(0).getModel().getFrequencyModel().getFrequency(0), 1e-6);
        assertEquals(1.004932107403357, result.get(1).getRate(), 1e-6);
    }
}
