package org.fhcrc.matsen.startree.gson

import org.scalatest.FlatSpec
import org.scalatest.matchers.ShouldMatchers
import org.fhcrc.matsen.startree.MatrixUtils._
import org.apache.commons.math.linear.ArrayRealVector
import org.fhcrc.matsen.startree.StarTreeTraces

/**
 * Created by cmccoy on 1/6/14.
 */
class StarTreeTracesTestCase extends FlatSpec with ShouldMatchers {
  "round trips" should "succeed" in {
    val r = 4
    val c = 12

    val cn1 = createRealMatrix(r, c, 1.0)
    val un1 = createRealMatrix(r, c, 12.1)
    val cs1 = createRealMatrix(r, c, 0.7)
    val us1 = createRealMatrix(r, c, 7.7)
    val bl1 = createRealMatrix(r, c, 14.1)

    for(i <- 0 until r) {
      bl1.setEntry(i, 0, 0)
      bl1.setEntry(i, 0, 1)
    }

    val state = new ArrayRealVector(Array[Double](1, 20, 40, 80));
    val cov1 = createRealVector(c, 1.0);

    val t1 = new StarTreeTraces(state, cn1, cs1, un1, us1, cov1.toArray(), bl1)

    // Serialize
    val gson = org.fhcrc.matsen.startree.gson.getGsonBuilder.create()

    val serString = gson.toJson(t1)

    val fromStr = gson.fromJson(serString, classOf[StarTreeTraces])

    cov1 getData() zip (fromStr getCoverage() getData) foreach {
      case (orig, d) =>
        orig should be (d plusOrMinus 1e-6)
    }

    fromStr.getConditionalNonsynonymous.getNorm should be (cn1.getNorm plusOrMinus 1e-6)
    fromStr.getUnconditionalNonsynonymous.getNorm should be (un1.getNorm plusOrMinus 1e-6)
    fromStr.getConditionalSynonymous.getNorm should be (cs1.getNorm plusOrMinus 1e-6)
    fromStr.getUnconditionalSynonymous.getNorm should be (us1.getNorm plusOrMinus 1e-6)
  }
}
