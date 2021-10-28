#!/bin/bash 
PATH_pdf=$1 # INPUT
PATH_png=$2

pdftoppm  ${PATH_pdf} ${PATH_png} -png
mv "${PATH_png}-1.png" "${PATH_png//-1.png/}"