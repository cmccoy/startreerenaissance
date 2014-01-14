package org.fhcrc.matsen.startree

import java.io._
import com.google.common.io.Files
import java.util.zip.{GZIPOutputStream, GZIPInputStream}

case class UpdateOldConfig(jsonPaths: Seq[File] = Seq())

/**
 * Created by cmccoy on 1/9/14.
 */
object UpdateOldJsonResults {
  val appName = "StarTreeUpdateOldJsonResults"

  private def transformName(path: String) : String = {
    val newPath = Files.getFileExtension(path) match {
      case "gz" => path.replaceAll("\\.json\\.gz", "_updated.json.gz")
      case "json" => path.replaceAll("\\.json\\.gz", "_updated.json.gz")
      case _ => throw new IllegalArgumentException("Unknown extension on: " + path)
    }

    require(newPath != path, "path unchanged!")

    newPath
  }

  private def inputStreamOfFile(file: File) : java.io.InputStream = {
    val buf = new BufferedInputStream(new FileInputStream(file))

   Files.getFileExtension(file.getAbsolutePath) match {
      case "gz" => new GZIPInputStream(buf)
      case _ => buf
    }
  }

  private def outputStreamOfFile(file: File) : java.io.OutputStream = {
    val buf = new BufferedOutputStream(new FileOutputStream(file))

    Files.getFileExtension(file.getAbsolutePath) match {
      case "gz" => new GZIPOutputStream(buf)
      case _ => buf
    }
  }

  val parser = new scopt.OptionParser[UpdateOldConfig](appName) {
    head(appName, "0.1")
    arg[File]("<json>...") unbounded() action { (x, c) => c.copy(jsonPaths=c.jsonPaths :+ x) } validate {
      x => if(x.isFile()) success else failure(x.toString() + " is not a file.")
    }
    help("help") text("prints this usage text")
  }

  def run(config: UpdateOldConfig) {
    config.jsonPaths foreach {
      jsonPath => {
        println(jsonPath)
        val inputStream = inputStreamOfFile(jsonPath)

        val reader = new InputStreamReader(inputStream)
        val deser = new com.google.gson.JsonParser()

        val g = gson.getGsonBuilder.create()
        val unsmoothed = g.fromJson(deser.parse(reader).getAsJsonObject.get("unsmoothed"), classOf[StarTreeTraces])
        unsmoothed.getSmoothed(true)

        reader.close()

        val r = new StarTreeRenaissanceResult("unknown", unsmoothed, true, 2000)

        if(Files.getFileExtension(jsonPath.getCanonicalPath) == "gz")
          require(jsonPath.getCanonicalPath.matches(".*.json.gz"))
        if(Files.getFileExtension(jsonPath.getCanonicalPath) == "json")
          require(jsonPath.getCanonicalPath.matches(".*.json"))

        val outputPath = transformName(jsonPath.getCanonicalPath)
        println("Writing to " + outputPath)

        val writer = new OutputStreamWriter(outputStreamOfFile(new java.io.File(outputPath)))
        g.toJson(r, writer)
        writer.close()
      }
    }
  }

  def main(args: Array[String]) {
    parser.parse(args, UpdateOldConfig()) map { config =>
      run(config)
    } getOrElse {
      // arguments are bad, error message will have been displayed
      System.exit(1)
    }
  }
}
