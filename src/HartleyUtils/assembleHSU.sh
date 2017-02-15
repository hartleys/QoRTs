#!/bin/bash

FULLDATE=$(date)
FULLDATESEC=$(date --date="$FULLDATE" +%s)

QORTS_DATE_LINE="  final val QORTS_COMPILE_DATE = \"$FULLDATE\"; \\/\\/ REPLACE_THIS_QORTS_DATE_VARIABLE_WITH_DATE          (note this exact text is used in a search-and-replace. Do not change it.)"
QORTS_TIME_LINE="  final val QORTS_COMPILE_TIME : Long = $FULLDATESEC; \\/\\/ REPLACE_THIS_QORTS_DATE_VARIABLE_WITH_TIME          (note this exact text is used in a search-and-replace. Do not change it.)"

sed -i -e "s/^.*REPLACE_THIS_QORTS_DATE_VARIABLE_WITH_DATE.*/$QORTS_DATE_LINE/g" ./src/main/scala/runner/runner.scala
sed -i -e "s/^.*REPLACE_THIS_QORTS_DATE_VARIABLE_WITH_TIME.*/$QORTS_TIME_LINE/g" ./src/main/scala/runner/runner.scala

if [ $# -gt 0 ]
then
  echo "Setting version $1"
  QORTS_VERS_LINE="  final val QORTS_VERSION = \"$1\"; \\/\\/ REPLACE_THIS_QORTS_VERSION_VARIABLE_WITH_VERSION_NUMBER          (note this exact text is used in a search-and-replace. Do not change it.)"
  sed -i -e "s/^.*REPLACE_THIS_QORTS_VERSION_VARIABLE_WITH_VERSION_NUMBER.*/$QORTS_VERS_LINE/g" ./src/main/scala/runner/runner.scala
fi

sbt < sbtAssemblyCommand.txt
