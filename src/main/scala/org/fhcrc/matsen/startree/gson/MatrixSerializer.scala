package org.fhcrc.matsen.startree.gson

import com.google.gson._
import java.lang.reflect.Type
import org.apache.commons.math.linear.{ArrayRealVector,BlockRealMatrix,RealMatrix,RealVector}

class RealMatrixSerializer extends JsonSerializer[RealMatrix] {
  def serialize(src: RealMatrix, typeOfSrc : Type, context: JsonSerializationContext) = {
    context.serialize(src.getData)
  }
}
object RealMatrixSerializer {
  val serializedType = classOf[RealMatrix]
}

class RealMatrixDeserializer extends JsonDeserializer[RealMatrix] {
  def deserialize(json: JsonElement, typeOfElement : Type, context: JsonDeserializationContext) = {
    val arrayValue : Array[Array[Double]] = context.deserialize(json, classOf[Array[Array[Double]]])
    new BlockRealMatrix(arrayValue)
  }
}
object RealMatrixDeserializer {
  val serializedType = classOf[RealMatrix]
}

class RealVectorSerializer extends JsonSerializer[RealVector] {
  def serialize(src: RealVector, typeOfSrc : Type, context: JsonSerializationContext) = {
    context.serialize(src.getData)
  }
}
object RealVectorSerializer {
  val serializedType = classOf[RealVector]
}

class RealVectorDeserializer extends JsonDeserializer[RealVector] {
  def deserialize(json: JsonElement, typeOfElement : Type, context: JsonDeserializationContext) = {
    val arrayValue : Array[Double] = context.deserialize(json, classOf[Array[Double]])
    new ArrayRealVector(arrayValue)
  }
}

object RealVectorDeserializer {
  val serializedType = classOf[RealVector]
}
