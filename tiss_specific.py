#!/usr/bin/python
#-*-coding : utf-8-*-
#Copyright(c) 2019 - leimingJiang <leiming8886@163.com>
import optparse
import pandas as pd
import re
import copy
#FileName = raw_input('Please enter your a file name: ')
parse=optparse.OptionParser()
parse.add_option('-f','--file',dest='file',action='store',metavar='input files',help='enter your transcript (contain the transcript ID)')
parse.add_option('-o','--out',dest='outfile',action='store',metavar='output files',help='assign your output file')
parse.add_option('-t','--tss',dest='tss',action='store',metavar='tissue name',help='please enter tissue name')
parse.add_option('-i','--t',dest='tiss',action='store_true',metavar='tissue name',help='please enter tissue name')
#store_true
(options,args) = parse.parse_args()
inPutFileName = options.file
outPutFileName = options.outfile
tss_name = options.tss
tss_Type = options.tiss
df = pd.read_csv("expression.txt",sep='\t')
col_name_del_probes_De=list(df.columns)
col_name_del_probes_De.remove("Probes")
col_name_del_probes_De.remove("Description")
'''
tiss=[ 'adipose', 'adrenal', 'brain', 'breast',
       'colon', 'heart', 'kidney', 'liver', 'lung', 'lymphNode', 'ovary',
       'prostate', 'skeltalMuscle', 'whiteBloodCell', 'testes', 'thyroid',
       'testes_R', 'brain_R', 'placenta_R', 'foreskin_R', 'hLF_r2', 'hLF_r1']

'''
dict_tiss={}

# print(list(df.loc[df['adipose'] != 0,'Probes']))
#dict_tiss dict key tiss,value probes
for col_line in col_name_del_probes_De:
    dict_tiss[col_line] = list(df.loc[df[col_line] != 0,'Probes'])
#test 'adipose'
# pick='adipose'
pick=tss_name
Pick_pro=dict_tiss[pick]
pick_sp=copy.deepcopy(dict_tiss)
del pick_sp[pick]
lncRNA_left_pro=set()
for col_line in pick_sp:
    for probe in pick_sp[col_line]:
        lncRNA_left_pro.add(probe)
# print(len(list(lncRNA_left_pro)))
# print ('XLOC_000397'in lncRNA_left_pro)
# print('XLOC_000397'in Pick_pro)
lncRNA_left_pro.remove('XLOC_000397')
Pick_pro.remove('XLOC_000397')
print (set(Pick_pro)-lncRNA_left_pro)
#hg19tohg38 del XLOC_000397
#set_lncRNA key probe,value temp_list[chr,start,end,strand]
ref_gtf=open('CrossMap_CabiliSuppDataSet1_lincRNAs_transcripts.gtf','r')
fr=ref_gtf.readlines()
set_lncRNA={}
for line in fr:
    # key=gene_ID,value  ruturn chr,start,end,strand
    line=line.strip()
    arr_line=line.split('\t')
    chr=arr_line[0]
    start=int(arr_line[3])
    end=int(arr_line[4])
    strand=arr_line[6]
    temp_list=[chr,start,end,strand]
    last_arr_tem=arr_line[-1].split(';')[0]
    last_arr_tem=last_arr_tem.replace('\"','')
    last_arr_id=last_arr_tem.split(' ')[1]
    # print(last_arr_id)
    if last_arr_id not in set_lncRNA:
        set_lncRNA[last_arr_id]=[]
        set_lncRNA[last_arr_id].append(temp_list)
    else:
        if temp_list in set_lncRNA[last_arr_id]:
            continue
        else:
            set_lncRNA[last_arr_id].append(temp_list)
ref_gtf.close()
# print(set_lncRNA["XLOC_013608"])


#open gtf to align
#the first five line of gtf_union.gtf is test,3 overlop,2 return
sam_gtf=open(inPutFileName,'r')
sam_gtf_out=open(outPutFileName,'w')
# def return_line(*list):
#     chr,start,end, strand=list

for line in sam_gtf.readlines():
    line=line.strip()
    arr_line=line.split('\t')
    chr='chr'+arr_line[0]
    start=int(arr_line[3])
    end=int(arr_line[4])
    strand=arr_line[6]
    #temp_list=[chr,start,end,strand]
    num_temp_list=0
    for key_id in set_lncRNA:
        for list in set_lncRNA[key_id]:
            if list[-1]==strand and list[0] == chr:
                if (start- list[1])*(list[2]-start)>=0 or (end- list[1])*(list[2]-end)>=0:
                    num_temp_list+=1
    if num_temp_list>0:
        continue
    else:
        sam_gtf_out.write(line+'\n')
sam_gtf.close()
sam_gtf_out.close()


# #open gtf to align
# #the first five line of gtf_union.gtf is test,3 overlop,2 return
sam_gtf=open(inPutFileName,'r')
sam_gtf_out=open(outPutFileName,'w')
# def return_line(*list):
#     chr,start,end, strand=list

for line in sam_gtf.readlines():
    line=line.strip()
    arr_line=line.split('\t')
    chr='chr'+arr_line[0]
    start=int(arr_line[3])
    end=int(arr_line[4])
    strand=arr_line[6]
    #temp_list=[chr,start,end,strand]
    num_temp_list=0
    num_list_pick=0
    for key_id in lncRNA_left_pro:
        for list in set_lncRNA[key_id]:
            if list[-1]==strand and list[0] == chr:
                if (start- list[1])*(list[2]-start)>=0 or (end- list[1])*(list[2]-end)>=0:
                    num_temp_list+=1
    for key_id in Pick_pro:
        for list in set_lncRNA[key_id]:
            if list[-1]==strand and list[0] == chr:
                if (start- list[1])*(list[2]-start)>=0 or (end- list[1])*(list[2]-end)>=0:
                    num_list_pick+=1
    if num_temp_list>0:
        continue
    else:
        if num_list_pick>0:
            sam_gtf_out.write(line+'\n')
        else:
            continue
sam_gtf.close()
sam_gtf_out.close()




