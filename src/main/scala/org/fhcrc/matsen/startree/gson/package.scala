package org.fhcrc.matsen.startree

import com.google.gson.GsonBuilder

/**
 * Created by cmccoy on 1/6/14.
 */
package object gson {
  def getGsonBuilder : GsonBuilder = {
    new GsonBuilder()
      .registerTypeAdapter(RealMatrixSerializer.serializedType, new RealMatrixSerializer)
      .registerTypeAdapter(RealMatrixDeserializer.serializedType, new RealMatrixDeserializer)
      .registerTypeAdapter(RealVectorSerializer.serializedType, new RealVectorSerializer)
      .registerTypeAdapter(RealVectorDeserializer.serializedType, new RealVectorDeserializer)
      .serializeSpecialFloatingPointValues()
  }
}
