
#!/usr/bin/python
#-*-coding : utf-8-*-
#Copyright(c) 2019 - leimingJiang <leiming8886@163.com>
import optparse
import re
######################### define input and output######################################
#FileName = raw_input('Please enter your a file name: ')
parse=optparse.OptionParser()
parse.add_option('-f','--file',dest='file',action='store',metavar='input files',help='enter your transcript (contain the transcript ID)')
parse.add_option('-o','--out',dest='outfile',action='store',metavar='output files',help='assign your output file')
parse.add_option('-g','--gtf',dest='gtf',action='store',metavar='gtf file name',help='please enter your gtf files')

(options,args) = parse.parse_args()
inPutFileName = options.file
outPutFileName = options.outfile
FileType = options.gtf
# inPutFileName = "intersect_plek_cnci.txt"
# outPutFileName = "k510_ouy.gtf"
# FileType =  "K510.gtf"
#set lncRNA transcript ID;
fr_lnc= open(inPutFileName)
#fr_lnc_firs=fr_lnc.readline()
lncRNA_set={}
for line in fr_lnc.readlines():
	line.strip()
	arr=line.split('\t')
	if arr[0] == "transcript ID":
		continue
	lncRNA_set[arr[0]]=arr[1]
fr_lnc.close()
#print([key for key in lncRNA_set])
#open GTF
fr_gtf = open(FileType)
fr_out = open(outPutFileName,'w')
for line in fr_gtf.readlines():
	line.strip()
	arr=line.split('\t')
	array=arr[-1].split(';')
	del array[-1]
	for line1 in array:
		line1=line1.strip()
		line1=line1.replace('\"','')
		#print(line1)
		line1 = line1.strip()
    # \s is a special regular expression character class for white spaces. str.split does not accept regular expressions. If you want to split with regular expressions, you have to use re.split.
		arr1=line1.split(' ')
		if arr1[0]=="transcript_id" and lncRNA_set.get(arr1[1]):
			#print(line1)
			fr_out.write(line)
fr_gtf.close()
fr_out.close()
