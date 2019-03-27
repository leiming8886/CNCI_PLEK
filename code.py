#!/usr/bin/python
#-*-coding : utf-8-*-
#Copyright(c) 2019 - leimingJiang <leiming8886@163.com>
from Bio import SeqIO
#import pandas as pd
import optparse
import os
######################### define input and output######################################
#FileName = raw_input('Please enter your a file name: ')
parse=optparse.OptionParser()
parse.add_option('-f','--file',dest='file',action='store',metavar='input files',help='enter your transcript (sequence or gtf)')
#parse.add_option('-o','--out',dest='outfile',action='store',metavar='output files',help='assign your output file')
parse.add_option('-p','--parallel',dest='parallel',action='store',metavar='prallel numbers',help='please enter your specified speed ratio')
parse.add_option('-m','--model',dest='model',action='store',metavar='model types',default='ve',help='please enter your specified classification model')
parse.add_option('-g','--gtf',dest='gtf',action='store_true',metavar='gtf file name',help='please enter your gtf files')
parse.add_option('-d','--directory',dest='directory',action='store',metavar='',help='if your input file is gtf type please enter RefGenome directory')
parse.add_option('-i','--cnci',dest='cnci',action='store',metavar='',help='enter path of the CNCI.py')
parse.add_option('-k','--plek',dest='plek',action='store',metavar='',help='enter path of the PLEK.py')

(options,args) = parse.parse_args()
inPutFileName = options.file
#outPutFileName = options.outfile
Parallel = options.parallel
ClassModel = options.model
FileType = options.gtf
Directory = options.directory
CNCI = options.cnci #set absolute path
PLEK = options.plek #set absolute path

CNCIPATH=os.path.split(os.path.realpath(CNCI))[0]
PLEK=os.path.split(os.path.realpath(PLEK))[0]
outPutFileName=os.path.splitext(inPutFileName)[0]
def sub_array(A,B):
    x=set(A)
    y=set(B)
    return list(x - y)

def intersect_array(A,B):
    x=set(A)
    y=set(B)
    return list(x & y)


def union_array(A,B):
    x=set(A)
    y=set(B)
    return list(x | y)

def TwoLineFasta (Seq_Array):
    Tmp_sequence_Arr = []
    Tmp_trans_str = ''
    for i in range(len(Seq_Array)):
        Seq_Array[i]=Seq_Array[i].strip()
        if '>' in Seq_Array[i]:
            if i == 0:
                Tmp_sequence_Arr.append(Seq_Array[i])
            else:
                Tmp_sequence_Arr.append(Tmp_trans_str)
                Tmp_sequence_Arr.append(Seq_Array[i])
                Tmp_trans_str = ''
        else:
            if i == len(Seq_Array) - 1:
                Tmp_trans_str = Tmp_trans_str + str(Seq_Array[i])
                Tmp_sequence_Arr.append(Tmp_trans_str)
            else:
                Tmp_trans_str = Tmp_trans_str + str(Seq_Array[i])
    return Tmp_sequence_Arr

def noTwolineFasta (Seq_Array,wid=50):
	Tmp_sequence_Arr = []
	width=int(wid)
	for i in range(len(Seq_Array)):
		Seq_Array[i]=Seq_Array[i].strip()
		len_fa=len(Seq_Array[i])
		start = 0
		end = start + width
		if '>' in Seq_Array[i]:
			Tmp_sequence_Arr.append(Seq_Array[i])
		else:
			while end < len_fa:
			#output=substr(seq_record.seq, i, 50)
				Tmp_sequence_Arr.append(str(Seq_Array[i][start:end]))
				start = start + width
				end += width
			if end >= len_fa:
				Tmp_sequence_Arr.append(str(Seq_Array[i][start:]))
	return Tmp_sequence_Arr

'''file1s=`ls *.gtf`
for file1 in ${file1s};
do python /home/lmjiang/software/CNCI/CNCI.py -f $file1 -g -o ${file1%.gtf} -m ve -p 8 -d /home/lmjiang/software/CNCI/hg38.2bit ;done
'''
#CNCI and PLEK code
if FileType:
	os.system('python ' + CNCIPATH + '/CNCI.py -f '+inPutFileName+' -g -o '+outPutFileName+' -m ve -p '+Parallel+' -d ' +Directory)
#fasta is not TwoLineFasta, fastaFiles = inPutFileName + '.fa', so need to convert format TwoLineFasta
	fastaFiles = outPutFileName + '.gtf.fa'
	fastaFiles_twoline=outPutFileName +'_plek'+'.fa'
	GtfInFiles = open(fastaFiles)
	inFilesArr = GtfInFiles.read()
	sequence_Arr = inFilesArr.split('\n')
	sLen = len(sequence_Arr) - 1#the last row is the null due to split (\n)
	del sequence_Arr[sLen]
	ARRAY =  TwoLineFasta(sequence_Arr)
	fr = open(fastaFiles_twoline,'w')
	for line in ARRAY:
		fr.write(line+"\n")
	fr.close()
	os.system('python ' + PLEK + '/PLEK.py '+' -fasta '+fastaFiles_twoline+' -out '+outPutFileName+'_PLEK'+' -thread 10 ')
else:
	os.system('python ' + PLEK + '/PLEK.py '+' -fasta '+inPutFileName+' -out '+outPutFileName+'_PLEK'+' -thread 10 ')
	fastaFiles_notwoline=outPutFileName+"_notwoline"+'.fa'
	GtfInFiles = open(inPutFileName)
	inFilesArr = GtfInFiles.read()
	sequence_Arr = inFilesArr.split('\n')
	sLen = len(sequence_Arr) - 1#the last row is the null due to split (\n)
	del sequence_Arr[sLen]
	ARRAY =  noTwolineFasta(sequence_Arr,50)
	fr = open(fastaFiles_notwoline,'w')
	for line in ARRAY:
		fr.write(line+"\n")
	fr.close()
	os.system('python ' + CNCIPATH + '/CNCI.py -f '+fastaFiles_notwoline+' -o '+outPutFileName+' -m ve -p '+Parallel)

#set hash of two output of the out_plek and out_cnci
#[key for key in d]
out_plek=outPutFileName+'_PLEK'
out_cnci=outPutFileName+'/CNCI.index'
#out_plek = "K510"
#out_cnci = "CNCI.index"
#set hash of two output of the out_plek and out_cnci
#[key for key in d]
out_set_plek={}
out_set_cnci={}
pl_fr = open(out_plek)
cn_fr = open(out_cnci)
#g_pl_fr1 = pl_fr.readline()
for line in  pl_fr.readlines():
    line=line.replace('>','')
    line1=line.strip()
    line_c=line1.split('\t')
    #print(line_c[2])
    if line_c[0] == 'Coding':
        continue
    out_set_plek[line_c[2]]=line_c[0]+"\t"+line_c[1]
#d.items()
gi_p=[key for key in out_set_plek]
#print(gi_p)

for line in  cn_fr.readlines():
    line1=line.strip()
    line_c=line1.split('\t')
    if line_c[1] ==  'coding' or line_c[1] ==  'index':
        continue
    out_set_cnci[line_c[0]]=line_c[1]+"\t"+line_c[2]
#d.items()
gi_c=[key for key in out_set_cnci]
#print(gi_c)
pl_fr.close()
cn_fr.close()
union=open(outPutFileName+"union_plek_cnci.txt",'w')
inter=open(outPutFileName+"intersect_plek_cnci.txt",'w')
union.write('transcript ID\tPLEK_index\tPLEK_score\tCNCI_index\tCNCI_score\n')
inter.write('transcript ID\tPLEK_index\tPLEK_score\tCNCI_index\tCNCI_score\n')
for key in union_array(gi_p,gi_c):
	if not out_set_cnci.get(key):
		out_set_cnci[key]="null\tnull"
	if not out_set_plek.get(key):
		out_set_plek[key]="null\tnull"
	union.write(key+'\t'+out_set_plek[key]+'\t'+out_set_cnci[key]+'\n')
	#print(key+'\tnonconding\n')
for key in intersect_array(gi_p,gi_c):
	inter.write(key+'\t'+out_set_plek[key]+'\t'+out_set_cnci[key]+'\n')
	#print(key + '\tnonconding\n')
union.close()
inter.close()
#run gtf pass : python3 PLEK_CNCI.py -f K510.gtf -p 6 -m ve -g -d /home/lmjiang/software/CNCI/hg38.2bit -i /home/lmjiang/software/CNCI/CNCI.py -k /home/lmjiang/software/PLEK.1.2/PLEK.py
#run fasta pass: python3 PLEK_CNCI.py -f K510.fa -p 6 -m ve -i /home/lmjiang/software/CNCI/CNCI.py -k /home/lmjiang/software/PLEK.1.2/PLEK.py
