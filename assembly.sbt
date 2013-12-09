import AssemblyKeys._ // put this at the top of the file

assemblySettings

mainClass := Some("org.fhcrc.matsen.startree.spark.StarTreeSpark")

mergeStrategy in assembly <<= (mergeStrategy in assembly) { (old) =>
  {
    case PathList("dr", "math", xs @ _*) => MergeStrategy.first
    case PathList("dr", xs @ _*) => MergeStrategy.last
    case x => old(x)
  }
}
