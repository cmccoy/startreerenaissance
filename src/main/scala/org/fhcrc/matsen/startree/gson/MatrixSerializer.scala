package org.fhcrc.matsen.startree.gson

import com.google.gson._
import java.lang.reflect.Type
import org.apache.commons.math.linear.{ArrayRealVector,BlockRealMatrix,RealMatrix,RealVector};
import org.fhcrc.matsen.startree._

class BlockRealMatrixSerializer extends JsonSerializer[BlockRealMatrix] {
  def serialize(src: BlockRealMatrix, typeOfSrc : Type, context: JsonSerializationContext) = {
    context.serialize(src.getData)
  }
}

object BlockRealMatrixSerializer {
  val serializedType = classOf[RealMatrix];
}

class BlockRealMatrixDeserializer extends JsonDeserializer[RealMatrix] {
  def deserialize(json: JsonElement, typeOfElement : Type, context: JsonDeserializationContext) = {
    val arrayValue : Array[Array[Double]] = context.deserialize(json, classOf[Array[Array[Double]]])
    new BlockRealMatrix(arrayValue)
  }
}

class ArrayRealVectorDeserializer extends JsonDeserializer[RealVector] {
  def deserialize(json: JsonElement, typeOfElement : Type, context: JsonDeserializationContext) = {
    val arrayValue : Array[Double] = context.deserialize(json.getAsJsonObject().get("data"), classOf[Array[Double]])
    new ArrayRealVector(arrayValue)
  }
}

object ArrayRealVectorDeserializer {
  val serializedType = classOf[RealVector]
}
