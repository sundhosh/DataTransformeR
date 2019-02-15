Publish
R package Publish

Installation
To install the development version of Publish run the following commands from within R

library(devtools)
install_github('sundhosh/DataTransformeR')
Trouble shooting
To install a package from github you need a program to unzip the download. If you don't have such a program and the install_github above command failed, then you should try

library(devtools)
options(unzip="internal")
install_github('sundhosh/DataTransformeR')
