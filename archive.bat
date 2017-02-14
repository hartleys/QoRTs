@ECHO OFF

SET /P VER=Please enter the OLD version number: 

echo "COPYING OLD FILES..."
:: go to the directory containing this .bat file:
cd %~dp0

mkdir archive\v%VER%
xcopy /E /I /Q doc archive\v%VER%\doc
xcopy /E /I /Q jarHtml archive\v%VER%\jarHtml
xcopy /E /I /Q jarMd archive\v%VER%\jarMd
xcopy /E /I /Q Rhtml archive\v%VER%\Rhtml
xcopy /E /I /Q stylesheets archive\v%VER%\stylesheets

xcopy * archive\v%VER%\

SET /P VER=Script complete. Hit enter to end.
