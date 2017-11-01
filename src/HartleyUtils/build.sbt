import AssemblyKeys._

assemblySettings

scalaVersion := "2.12.4"

scalacOptions := Seq("-unchecked", "-deprecation")

jarName in assembly := "QoRTs.jar"

mainClass in assembly := Some("runner.runner")

