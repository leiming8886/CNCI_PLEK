#!/bin/sh
#python 
tar zvxf PLEK.1.2.tar.gz 
cd PLEK.1.2
python PLEK_setup.py
cd ../
upzip CNCI-master.zip
cd CNCI-master
unzip libsvm-3.0.zip
cd libsvm-3.0
make
cd ../../
#source setup.sh
