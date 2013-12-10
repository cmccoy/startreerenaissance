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
    "com.novocode" % "junit-interface" % "0.9" % "test",
    "org.scalatest" %% "scalatest" % "1.9.1" % "test",
    "com.github.scopt" %% "scopt" % "3.2.0",
    "org.apache.spark" % "spark-core_2.9.3" % "0.8.0-incubating" % "provided"
  )
}

runMain in Compile <<= Defaults.runMainTask(fullClasspath in Compile, runner in (Compile, run))

run in Compile <<= Defaults.runTask(fullClasspath in Compile, mainClass in (Compile, run), runner in (Compile, run)) 

