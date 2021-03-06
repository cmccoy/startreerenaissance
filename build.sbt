name := "startreerenaissance"

version := "0.1"

scalaVersion := "2.9.3"

net.virtualvoid.sbt.graph.Plugin.graphSettings

resolvers ++= Seq(
  "Akka Repository" at "http://repo.akka.io/releases/",
  "Spray Repository" at "http://repo.spray.cc/",
  "snapshots" at "http://oss.sonatype.org/content/repositories/snapshots",
  "releases" at "http://oss.sonatype.org/content/repositories/releases")



libraryDependencies ++= {
  Seq(
    "com.novocode" % "junit-interface" % "0.10" % "test",
    "org.scalatest" %% "scalatest" % "1.9.1" % "test",
    "org.apache.commons" % "commons-math3" % "3.2",
    "com.github.scopt" %% "scopt" % "3.2.0",
    "com.amazonaws" % "aws-java-sdk" % "1.4.7",
    "org.apache.spark" % "spark-core_2.9.3" % "0.8.1-incubating" % "provided"
  )
}

runMain in Compile <<= Defaults.runMainTask(fullClasspath in Compile, runner in (Compile, run))

run in Compile <<= Defaults.runTask(fullClasspath in Compile, mainClass in (Compile, run), runner in (Compile, run))

javacOptions in (Compile, compile) += "-Xlint"

fork in Test := true

addCommandAlias("runLocal", ";compile" +
                ";runMain org.fhcrc.matsen.startree.StarTreeMain simulate/test.json simulate/gy94_mixture.fasta simulate/gy94_mixture.bam gy94_mixture.log")
