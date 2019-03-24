#!/usr/bin/python
#-*-coding : utf-8-*-
#Copyright(c) 2019 - leimingJiang <leiming8886@163.com>
from Bio import SeqIO
import pandas as pd
######################### define input and output######################################
#FileName = raw_input('Please enter your a file name: ')
parse=optparse.OptionParser()
parse.add_option('-f','--file',dest='file',action='store',metavar='input files',help='enter your transcript (sequence or gtf)')
parse.add_option('-o','--out',dest='outfile',action='store',metavar='output files',help='assign your output file')
parse.add_option('-p','--parallel',dest='parallel',action='store',metavar='prallel numbers',help='please enter your specified speed ratio')
parse.add_option('-m','--model',dest='model',action='store',metavar='model types',default='ve',help='please enter your specified classification model')
parse.add_option('-g','--gtf',dest='gtf',action='store_true',metavar='gtf file name',help='please enter your gtf files')
parse.add_option('-d','--directory',dest='directory',action='store',metavar='',help='if your input file is gtf type please enter RefGenome directory')
parse.add_option('-i','--cnci',dest='cnci',action='store',metavar='',help='enter path of the CNCI.py')
parse.add_option('-k','--plek',dest='plek',action='store',metavar='',help='enter path of the CNCI.py')

(options,args) = parse.parse_args()
inPutFileName = options.file
outPutFileName = options.outfile
Parallel = options.parallel
ClassModel = options.model
FileType = options.gtf
Directory = options.directory
CNCI = options.cnci #set absolute path
PLEK = options.plek #set absolute path

CNCIPATH=os.path.split(os.path.realpath(CNCI))[0]
PLEK=os.path.split(os.path.realpath(PLEK))[0]

outPutFileName=os.path.splitext(inPutFileName)[0]#input是文件名，out是没有后缀的名称


def sub_array(A,B):
    x=Set(A)
    y=Set(B)
    return list(x - y)

def intersect_array(A,B):
    x=Set(A)
    y=Set(B)
    return list(x & y)


def union_array(A,B):
    x=Set(A)
    y=Set(B)
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

def noTwolineFasta (Seq_Array,wid):
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
	fastaFiles = outPutFileName + '.fa'
	fastaFiles_twoline=outPutFileName +'plek'+'.fa'
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
	os.system('python ' + PLEK + '/PLEK.py '+' -fasta '+fastaFiles_twoline+' -out '+outPutFileName+' -thread 10 ')
else:
    os.system('python ' + PLEK + '/PLEK.py '+' -fasta '+inPutFileName+' -out '+outPutFileName+' -thread 10 ')
	fastaFiles_notwoline=outPutFileName+"_notwoline"+'.fa'
	GtfInFiles = open(inPutFileName)
	inFilesArr = GtfInFiles.read()
	sequence_Arr = inFilesArr.split('\n')
	sLen = len(sequence_Arr) - 1#the last row is the null due to split (\n)
	del sequence_Arr[sLen]
	ARRAY =  noTwolineFasta(sequence_Arr)
	fr = open(fastaFiles_notwoline,'w')
	for line in ARRAY:
		fr.write(line+"\n")
	fr.close()
	os.system('python ' + CNCIPATH + '/CNCI.py -f '+fastaFiles_notwoline+' -o '+outPutFileName+' -m ve -p '+Parallel)

out_plek=outPutFileName
out_cnci=outPutFileName+'/CNCI.index'

#set hash of two output of the out_plek and out_cnci
#[key for key in d]
out_set_plek={}
out_set_cnci={}
pl_fr=pd.read_table(out_set_plek)
cn_fr=pd.read_table(out_set_cnci)
