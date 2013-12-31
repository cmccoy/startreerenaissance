package org.fhcrc.matsen.startree.spark;

import com.esotericsoftware.kryo.Kryo;
import dr.evolution.alignment.SimpleAlignment;
import org.fhcrc.matsen.startree.TwoTaxonResult;

/**
 * Created by cmccoy on 12/9/13.
 */
public class StarTreeKryoRegistrator implements org.apache.spark.serializer.KryoRegistrator {
    @Override
    public void registerClasses(Kryo kryo) {
        kryo.register(SimpleAlignment.class);
        kryo.register(org.apache.commons.math.linear.ArrayRealVector.class);
        kryo.register(org.apache.commons.math.linear.BlockRealMatrix.class);
        kryo.register(TwoTaxonResult.class);
    }
}
