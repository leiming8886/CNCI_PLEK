### Install CNCI_PLEK
At the first time to running CNCI and PLEK, we suggest you to install "libsvm-3.0" that stored in our package.

```
git clone https://github.com/leiming8886/CNCI_PLEK.git
cd CNCI_PLEK
source setup.sh
```
### HELP for CNCI_PLEK subroutines

**INCP.py: an integration classification tool of CNCI and PLEK for identify coding or non-coding transcripts (fasta file and gtf file)**

#### Usage: INCP.py -i input_gtf -p parallel -g -r ref_2bit

Parameters:

 -i or --input : input file of fasta file or gtf file, if the input is fasta file,the file format must be the twolineFasta

 -p or --parallel : assign the running CUP numbers

 -g or --gtf : if you input file is gtf format, please use this parameter

 -r or --reference : if you use the -g  this parameter must be assigned, within this parameter please assign the path of your reference genome. Some reference files which has been prepared could be download at [hg38](hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit), [hg19](hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit).

Output: mainly contains 5 files

 input_no_suffix directory: the name is the intput file without suffix. it contain the file of CNCI.index, which is the result of the software CNCI output.

 input_no_suffix_PLEK: which is the result of the software PLEK output.

 input_no_suffix_union_plek_cnci.txt: the output of union of the software CNCI and PLEK, in which the first column is the tanscript ID

 input_no_suffix_intersect_plek_cnci.txt: the output of intersect of the software CNCI and PLEK, in which the first column is the tanscript ID
 
 venny_plek_cnci.pdf : the summary of the venny between the output of the CNCI and PLEK
#### venny image
<div align="center"><img width="390" height="285" src="https://github.com/leiming8886/picture/blob/master/venny.png"></div>
<p align="center">
Venny
</p>

**extract_lncRNA_gtf.py: A tool that extract lncRNA information of GTF format based on the tanscript ID of the candidate lncRNA**

#### Usage: extract_lncRNA_gtf.py -i input -g GTF -o output_name

Parameters:

 -i or --input : input file of the candidate lncRNA, in which the first column is the tanscript ID of the candidate lncRNA. This file also can be the output file of INCP.py

 -g or --gtf : GTF file corresponding to fasta in the INCP.py, where the last column contain the tanscript ID


 -o or --out : output name extracted lncRNA information of GTF format


**tiss_specific.py : A tool that extract cancer-specific lncRNA information of GTF format**

#### Usage: tiss_specific.py -i input -t breast -r hg38 -o output_name 

Parameters:

 -i or --input : input files of the candidate lncRNA gtf format, if the input files have two splited by ',', the first set control sample, the other set cancer sample. If the input files have only one, there are two situations. The one is based on a background control tissue. The other have not a background control tissue. Related knowledge can refer to the https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3185964/

 -t or --tissue : set tissue name if input files have the only one. If there are two input files, this parameter can be ignored.The background control tissue : 'adipose', 'adrenal', 'brain', 'breast', 'colon', 'heart', 'kidney', 'liver', 'lung', 'lymphNode', 'ovary', 'prostate', 'skeltalMuscle', 'whiteBloodCell', 'testes', 'thyroid', 'placenta', 'foreskin', 'hLF'. 

 -r or --reference : set refgene name,'hg38' and 'hg19' can be chosen.

 -o or --out : output name extracted lncRNA-specific information of GTF format

**digitalTochrdigital.py : A tool that check if the first column chromosome is a numeric value of chr+digital  or of GTF format**

#### Usage: digitalTochrdigital.py -i inputdir -o out_dir

Parameters:

 -i, --inputdir : input files of the gtf format

 -o or --out : Outdir, create if the directory does not exist


## EXAMPLE
you can use CNCI_PLEK subroutines like our example:

```
python INCP.py -i candidate.gtf -p 6 -g -r hg38.2bit
or 
python INCP.py -i candidate.fasta -p 6
python extract_lncRNA_gtf.py -i input -g candidate.gtf -o outname
python tiss_specific.py -i control.gtf,sample.gtf -r hg38 -o outname
or 
python tiss_specific.py -i sample.gtf -t breast -r hg38-o outname
```
