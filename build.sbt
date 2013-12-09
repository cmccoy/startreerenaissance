name := "startreerenaissance"

version := "0.1"

scalaVersion := "2.9.3"

resolvers += "Typesafe Repository" at "http://repo.typesafe.com/typesafe/releases/"

resolvers += "Spray Repository" at "http://repo.spray.cc/"

libraryDependencies += "com.novocode" % "junit-interface" % "0.9" % "test"

libraryDependencies += "org.apache.spark" % "spark-core_2.9.3" % "0.8.0-incubating" % "provided"
