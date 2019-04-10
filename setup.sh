#!/bin/sh
#python 
tar zvxf PLEK.1.2.tar.gz 
cd PLEK.1.2
python PLEK_setup.py
cd ../
upzip fileName.zip
cd CNCI-master
chmod ugo+x  twoBitToFa
unzip libsvm-3.0.zip
cd libsvm-3.0
make
cd ../../
#source setup.sh
<<<<<<< HEAD
#matplotlib-venn
=======
>>>>>>> 7e0eb251e058e81aa098e1ab6463a87a151b216b
