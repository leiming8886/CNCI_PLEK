### Install CNCI_PLEK
At the first time to running CNCI and PLEK, we suggest you to install "libsvm-3.0" that stored in our package.

```
git clone https://github.com/leiming8886/CNCI_PLEK/
cd CNCI_PLEK
source setup.sh
```

## HELP for CNCI_PLEK subroutines

**main.py: an integration classification tool of CNCI and PLEK for identify coding or non-coding transcripts (fasta file and gtf file)**

#### Usage: main.py -f input_gtf -p parallel -g -d ref_2bit

Parameters:

 -f or --file : input file of fasta file or gtf file, if the input is fasta file,the file format must be the twolineFasta

 -p or --parallel : assign the running CUP numbers

 -g or --gtf : if you input file is gtf format please use this parameter

 -d or --directory : if you use the -g  this parameter must be assigned, within this parameter please assign the path of your reference genome.Some reference files which has been prepared could be download at [hg38](hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit), [hg19](hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit)

Output: mainly contains 4 files

 input_no_suffix directory: the name is the iutput file name without suffix. it contain the file of CNCI.index, which is the result of the software CNCI output.

 input_no_suffix_PLEK: which is the result of the software PLEK output.

 input_no_suffix_union_plek_cnci.txt: the output of union of the software CNCI and PLEK, in which the first column is the tanscript ID

 input_no_suffix_intersect_plek_cnci.txt: the output of intersect of the software CNCI and PLEK, in which the first column is the tanscript ID


** extract_lncRNA_gtf.py: A tool that extract lncRNA information of GTF format based on the tanscript ID of the candidate lncRNA**

#### Usage: extract_lncRNA_gtf.py -f input -g GTF -o out_dir

Parameters:

 -f or --file : input file of the candidate lncRNA, in which the first column is the tanscript ID of the candidate lncRNA. This file also can be the output file of CNCI.py

 -g or --gtf : GTF file corresponding to fasta in the main.py, where the last column contain the tanscript ID


 -o or --out : OUTfile extract lncRNA information of GTF format


#### tiss_spe.py : 

## EXAMPLE
you can use CNCI_PLEK subroutines like our example:

```
python main.py -f candidate.gtf -p 6 -g -d hg38.2bit
or 
python main.py -f candidate.fasta -p 6

python extract_lncRNA_gtf.py -f test.index -g unannotation.gtf -o out
```
