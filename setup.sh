#!/bin/sh
#python 
tar zvxf PLEK.1.2.tar.gz 
cd PLEK.1.2
python PLEK_setup.py
cd ../
unzip fileName.zip
cd CNCI-master
chmod ugo+x  twoBitToFa
unzip libsvm-3.0.zip
cd libsvm-3.0
make
cd ../../
#source setup.sh
#matplotlib-venn
