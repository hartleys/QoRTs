

scalaVersion := "2.12.4"

scalacOptions := Seq("-unchecked", "-deprecation", "-feature")

assemblyJarName in assembly := "QoRTs.jar"

mainClass in assembly := Some("runner.runner")

