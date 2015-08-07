source ~/setall.sh

FULLDATE=$(date)

QORTS_DATE_LINE="  final val QORTS_COMPILE_DATE = \"$FULLDATE\"; \\/\\/ REPLACE_THIS_QORTS_DATE_VARIABLE_WITH_DATE          (note this exact text is used in a search-and-replace. Do not change it.)"
sed -i -e "s/^.*REPLACE_THIS_QORTS_DATE_VARIABLE_WITH_DATE.*/$QORTS_DATE_LINE/g" ./src/main/scala/runner/runner.scala

sbt < sbtAssemblyCommand.txt
