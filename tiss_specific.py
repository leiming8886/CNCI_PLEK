#!/usr/bin/python
#-*-coding : utf-8-*-
#Copyright(c) 2019 - leimingJiang <leiming8886@163.com>
import optparse
import pandas as pd
import re
import copy
#FileName = raw_input('Please enter your a file name: ')
parse=optparse.OptionParser()
parse.add_option('-f','--files',dest='files',type='string',action='store',metavar='input files',help='enter your transcript (contain the transcript ID)')
#,nargs="*"
parse.add_option('-o','--out',dest='outfile',type='string',action='store',metavar='output files',help='assign your output file')
parse.add_option('-t','--tss',dest='tss',type='string',action='store',metavar='tissue name',help='please enter tissue name')
parse.add_option('-g','--g',dest='hg',type='string',action='store',metavar='ref name',help='please enter hg38 or hg19')
(options,args) = parse.parse_args()
outPutFileName = options.outfile
tss_name = options.tss
refgene = options.hg
#inPutFileNames=options.files
# print(options.files,"\n")
inPutFileNames = list(options.files.split(","))
# print(inPutFileNames,'\n')
df = pd.read_csv("expression.txt",sep='\t')
col_name_del_probes_De=list(df.columns)
col_name_del_probes_De.remove("Probes")
col_name_del_probes_De.remove("Description")

'''
tiss=[ 'adipose', 'adrenal', 'brain', 'breast',
       'colon', 'heart', 'kidney', 'liver', 'lung', 'lymphNode', 'ovary',
       'prostate', 'skeltalMuscle', 'whiteBloodCell', 'testes', 'thyroid',
       'testes_R', 'brain_R', 'placenta_R', 'foreskin_R', 'hLF_r2', 'hLF_r1']
_R: from the dataset Rinn lab
'''
#set dic key = tiss_name. value = transcript IDs
dict_tiss={}
# print(list(df.loc[df['adipose'] != 0,'Probes']))
#dict_tiss dict key tiss,value probes
for col_line in col_name_del_probes_De:
    col_line_Dedu=col_line.split('_')[0]
    dict_tiss[col_line_Dedu] = list(df.loc[df[col_line] != 0,'Probes'])


if refgene == 'hg38':
    ref_gtf=open('lncRNA_hg38.gtf','r')
    print("hg38\n")
if refgene == 'hg19':
    ref_gtf=open('lncRNA_hg19.gtf','r')
    print("hg19")
fr=ref_gtf.readlines()
#set_lncRNA key probe,value temp_list[[chr,start,end,strand],[chr,start,end,strand]]
set_lncRNA={}
for line in fr:
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

if len(inPutFileNames)==1:
    # print("inPutFileNames ","1\n")
    sam_gtf = open(inPutFileNames[0], 'r')
    sam_gtf_out = open(outPutFileName, 'w')
    for line in sam_gtf.readlines():
        # line=sam_gtf.readline()
        line=line.strip()
        arr_line=line.split('\t')
        #chr='chr'+arr_line[0]
        chr= arr_line[0]
        start=int(arr_line[3])
        end=int(arr_line[4])
        strand=arr_line[6]
        num_temp_list=0
        if tss_name not in dict_tiss:
            for key_id in set_lncRNA:
                for list in set_lncRNA[key_id]:
                    if list[-1]==strand and list[0] == chr:
                        if (start- list[1])*(list[2]-start)>=0 or (end- list[1])*(list[2]-end)>=0:
                            num_temp_list+=1
            if num_temp_list>0:
                continue
            else:
                sam_gtf_out.write(line+'\n')
        else:
            # print(tss_name+'\n')
            num_list_pick = 0
            pick = tss_name
            Pick_pro = dict_tiss[pick]
            pick_sp = copy.deepcopy(dict_tiss)
            del pick_sp[pick]
            lncRNA_left_pro = set()
            #lncRNA_left_pro:del the tiss name, 剩余probes
            for col_line in pick_sp:
                for probe in pick_sp[col_line]:
                    lncRNA_left_pro.add(probe)
            # print(len(list(lncRNA_left_pro)))
            # print ('XLOC_000397'in lncRNA_left_pro)
            # print('XLOC_000397'in Pick_pro)
            # lncRNA_left_pro.remove('XLOC_000397')
            # Pick_pro.remove('XLOC_000397')
            #print(set(Pick_pro) - lncRNA_left_pro)
            # try:
            for key_id in lncRNA_left_pro:
                if key_id not in set_lncRNA:
                    continue
                for list in set_lncRNA[key_id]:
                    if list[-1] == strand and list[0] == chr:
                        if (start - list[1]) * (list[2] - start) >= 0 or (end - list[1]) * (list[2] - end) >= 0:
                            num_temp_list += 1
            for key_id in Pick_pro:
                if key_id not in set_lncRNA:
                    continue
                for list in set_lncRNA[key_id]:
                    if list[-1] == strand and list[0] == chr:
                        if (start - list[1]) * (list[2] - start) >= 0 or (end - list[1]) * (list[2] - end) >= 0:
                            num_list_pick += 1
            if num_temp_list > 0:
                continue
            else:
                if num_list_pick > 0:
                    sam_gtf_out.write(line + '\n')
    sam_gtf.close()
    sam_gtf_out.close()
if len(inPutFileNames)==2:
    sam_gtf1 = open(inPutFileNames[0], 'r')
    sam_gtf2 = open(inPutFileNames[1], 'r')
    set_gtf2 = set()
    set_gtf1 = set()
    for line in sam_gtf1:
        # key=gene_ID,value  ruturn chr,start,end,strand
        line = line.strip()
        arr_line = line.split('\t')
        chr = arr_line[0]
        start = int(arr_line[3])
        end = int(arr_line[4])
        strand = arr_line[6]
        temp_list = [chr, start, end, strand]
        last_arr_tem = arr_line[-1].split(';')[0]
        last_arr_tem = last_arr_tem.replace('\"', '')
        last_arr_id = last_arr_tem.split(' ')[1]
        set_gtf1.add(last_arr_id)
    sam_gtf1.close()
    for line in sam_gtf2 :
        # key=gene_ID,value  ruturn chr,start,end,strand
        line = line.strip()
        arr_line = line.split('\t')
        chr = arr_line[0]
        start = int(arr_line[3])
        end = int(arr_line[4])
        strand = arr_line[6]
        temp_list = [chr, start, end, strand]
        last_arr_tem = arr_line[-1].split(';')[0]
        last_arr_tem = last_arr_tem.replace('\"', '')
        last_arr_id = last_arr_tem.split(' ')[1]
        set_gtf2.add(last_arr_id)
    sam_gtf2.close()
    sam = open(inPutFileNames[1], 'r')
    sam_gtf_out = open(outPutFileName, 'w')
    set_left=set_gtf2-set_gtf1
    for line in sam :
        # key=gene_ID,value  ruturn chr,start,end,strand
        line = line.strip()
        arr_line = line.split('\t')
        chr = arr_line[0]
        start = int(arr_line[3])
        end = int(arr_line[4])
        strand = arr_line[6]
        temp_list = [chr, start, end, strand]
        last_arr_tem = arr_line[-1].split(';')[0]
        last_arr_tem = last_arr_tem.replace('\"', '')
        last_arr_id = last_arr_tem.split(' ')[1]
        if last_arr_id in set_left:
            sam_gtf_out.write(line + '\n')
    sam.close()
    sam_gtf_out.close()



