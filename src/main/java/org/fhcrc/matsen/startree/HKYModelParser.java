package org.fhcrc.matsen.startree;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Parameter;

import java.io.Reader;
import java.util.ArrayList;
import java.util.List;

/**
 * Parses JSON models
 * Created by cmccoy on 12/5/13.
 */
public class HKYModelParser {
    public static List<HKYAndRate> substitutionModel(final Reader reader) {
        JsonParser parser = new JsonParser();
        JsonObject parsed = parser.parse(reader).getAsJsonObject();

        JsonArray partitions = parsed.getAsJsonArray("partitions");
        List<HKYAndRate> result = new ArrayList<HKYAndRate>(partitions.size());
        for(int i = 0; i < partitions.size(); i++) {
            final JsonObject p = partitions.get(i).getAsJsonObject();
            final double kappa = p.get("parameters").getAsJsonObject().get("HKY85.kappa").getAsDouble();
            final double[] pi = new double[4];
            final JsonArray piNode = p.get("pi").getAsJsonArray();
            for(int j = 0; j < piNode.size(); j++) {
                pi[j] = piNode.get(j).getAsDouble();
            }

            assert p.get("name").getAsString().equals("HKY85");

            final double rate = p.get("rate").getAsJsonObject().get("Constant.value").getAsDouble();
            result.add(new HKYAndRate(kappa, pi, rate));
        }
        return result;
    }

    public static class HKYAndRate implements java.io.Serializable {
        private final double kappa;
        private final double[] frequencies;
        private final double rate;

        public HKYAndRate(double kappa, double[] frequencies, double rate) {
            this.kappa = kappa;
            this.frequencies = frequencies;
            this.rate = rate;
        }

        public double getKappa() {
            return kappa;
        }

        public double[] getFrequencies() {
            return frequencies;
        }

        public double getRate() {
            return rate;
        }

        public HKY getModel() {
            return new HKY(kappa, new FrequencyModel(Nucleotides.INSTANCE, frequencies));
        }

        public SiteRateModel getSiteRateModel() {
            return new GammaSiteRateModel(String.format("rate"),
                    new Parameter.Default(rate),
                    null,
                    1,
                    null);
        }
    }
}
