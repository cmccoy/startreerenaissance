package org.fhcrc.matsen.startree.spark;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import com.esotericsoftware.kryo.Kryo;
import dr.evolution.alignment.SimpleAlignment;
import org.fhcrc.matsen.startree.TwoTaxonResult;

/**
 * Created by cmccoy on 12/9/13.
 */
public class StarTreeKryoRegistrator implements org.apache.spark.serializer.KryoRegistrator {
    @Override
    public void registerClasses(Kryo kryo) {
        kryo.register(DenseDoubleMatrix2D.class);
        kryo.register(cern.colt.matrix.impl.RCDoubleMatrix2D.class);
        kryo.register(SimpleAlignment.class);
        kryo.register(TwoTaxonResult.class);
    }
}
