import AssemblyKeys._

assemblySettings

scalaVersion := "2.11.1"


jarName in assembly := "QoRTs.jar"

mainClass in assembly := Some("runner.runner")

