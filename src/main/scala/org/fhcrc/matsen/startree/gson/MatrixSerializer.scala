package org.fhcrc.matsen.startree.gson

import com.google.gson._
import java.lang.reflect.Type
import org.apache.commons.math.linear.BlockRealMatrix
import org.fhcrc.matsen.startree._

class BlockRealMatrixSerializer extends JsonSerializer[BlockRealMatrix] {
  def serialize(src: BlockRealMatrix, typeOfSrc : Type, context: JsonSerializationContext) = {
    context.serialize(src.getData)
  }
}

object BlockRealMatrixSerializer {
  val serializedType = classOf[BlockRealMatrix];
}

//class BlockRealMatrixDeserializer extends JsonDeserializer[BlockRealMatrix] {
  //def deserialize(json: JsonElement, typeOfElement : Type, context: JsonDeserializationContext) = {
    //val arrayValue = context.deserialize(json, classOf[Array[Array[Double]]])
    //new BlockRealMatrix(arrayValue)
  //}
//}
