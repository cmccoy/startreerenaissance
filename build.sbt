name := "startreerenaissance"

version := "0.1"

scalaVersion := "2.9.3"

resolvers ++= Seq(("Typesafe Repository" at "http://repo.typesafe.com/typesafe/releases/"),
  ("Spray Repository" at "http://repo.spray.cc/"))

libraryDependencies ++= Seq(("com.novocode" % "junit-interface" % "0.9" % "test"),
  ("org.apache.spark" % "spark-core_2.9.3" % "0.8.0-incubating" % "provided"))
