package org.fhcrc.matsen.startree.spark

import com.esotericsoftware.kryo.Kryo
import dr.evolution.alignment.SimpleAlignment
import org.fhcrc.matsen.startree.TwoTaxonResult
import org.apache.commons.math.linear.{ArrayRealVector,BlockRealMatrix}
import org.apache.spark.serializer.KryoRegistrator

/**
 * Created by cmccoy on 1/2/14.
 */
class StarTreeKryoRegistrator extends KryoRegistrator {
  // Classes to register with Kryo
  private[this] val classes = Seq(classOf[SimpleAlignment],
    classOf[ArrayRealVector],
    classOf[BlockRealMatrix],
    classOf[TwoTaxonResult])

  def registerClasses(kryo : Kryo) = {
    classes.foreach(kryo.register)
  }
}
