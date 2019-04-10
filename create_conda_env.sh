#!/bin/bash

conda create -y -n py2_7 python=2.7
source activate py2_7

conda install -y -n py2_7 numpy matplotlib
pip install matplotlib_venn
pip install matplotlib



