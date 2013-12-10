import AssemblyKeys._ // put this at the top of the file

assemblySettings

mainClass := Some("org.fhcrc.matsen.startree.spark.StarTreeSpark")

mergeStrategy in assembly <<= (mergeStrategy in assembly) { (old) =>
  {
    case PathList("dr", "math", xs @ _*) => MergeStrategy.first
    case PathList("dr", xs @ _*) => MergeStrategy.last
    case PathList("about.html") => MergeStrategy.discard
    case PathList("org", "w3c", xs @ _*) => MergeStrategy.first
    case x => old(x)
  }
}

assembleArtifact in packageScala := false

test in assembly := {}
