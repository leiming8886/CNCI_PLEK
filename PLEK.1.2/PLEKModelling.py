﻿#!/usr/bin/env python
#######################################################
#                                                               
#  PLEK Modelling - build predictor of lncRNAs and mRNAs based on k-mer scheme  
#  Authors: Aimin Li, Junying Zhang                               
#  Contacts: LiAiminMail@gmail.com, jyzhang@mail.xidian.edu.cn  
#  Webcite: https://sourceforge.net/projects/plek/              
#  Version: 1.2
#  Updated on: June 26, 2014  
#                                                               
#######################################################


__all__ = ['']
import os, sys, traceback, getpass, time, re, subprocess, commands
from threading import Thread


class GridOption: # set or get input parameters
	args=""
	
	def __init__(self, options):
		dirname = os.path.dirname(__file__)
		self.svmtrain_pathname = os.path.join(dirname, './svm-train')
		self.svmpredict_pathname = os.path.join(dirname, './svm-predict')
		self.svmscale_pathname = os.path.join(dirname, './svm-scale')
		self.mRNA_file = ""
		self.lncRNA_file = ""
		self.log2c = "0,5,1" # 0,1,2,3,4,5
		self.log2g = "0,-5,-1" # 0,-1,-2,-3,-4,-5
		self.prefix = "plek_mdl_"
		self.rangefile = ""
		self.is_mRNAlncRNA_balanced = 0
		self.thread_count = 12
		self.modelfile=""
		self.kmer=5
		self.min_seq_length=200
		self.unkown=0
		self.is_recompile=0		
		self.isoutmsg=0
		self.isrmtempfile=1
		self.isadjustweight=0
		self.nfold=10
		self.script_dir="./"
		self.islogging=1
		self.x_fold=1 # down-sampling fold
		self.parse_opts(options)

	def parse_opts(self, options):
		args=options # save options to args    
		if type(options) == str:
			options = options.split()
		i = 0
		pass_through_opts = []
		# get python script file's name, args[0]
		# determine it is with dir or only file_name
		script_name=args[0]    

		if script_name.rfind('/')>=0:
			self.script_dir=script_name[:script_name.rfind('/')+1]
    
			self.svmtrain_pathname = self.script_dir + 'svm-train'
			self.svmpredict_pathname = self.script_dir + 'svm-predict'
			self.svmscale_pathname = self.script_dir + 'svm-scale'
			
		while i < len(options):
			if options[i] == '-mRNA': # a file to provide mRNA transcript sequences in fasta format.
				i = i + 1
				self.mRNA_file = options[i]
			elif options[i] == '-lncRNA': # a file to provide lncRNA transcript sequences in fasta format.
				i = i + 1
				self.lncRNA_file = options[i]   
			elif options[i] == '-thread': # The number of threads for running the PLEK Modelling program. The bigger this number is, the faster it runs.
				i = i + 1
				self.thread_count = int(options[i])
			elif options[i] == '-log2c': # The search range of svm C parameter. Format: from, to, by. i.e. 1,5,1 means log2c=1,2,3,4,5.
				i = i + 1
				self.log2c = options[i]   
			elif options[i] == '-log2g': # The search range of svm G parameter. Format: from, to, by. i.e. 1,-4,-1 means log2g=1,0,-1,-2,-3,-4.
				i = i + 1
				self.log2g = options[i]
			elif options[i] == '-minlength': # The minimum length of sequences. The sequences whose lengths are more than minlength will be processed.
				i = i + 1
				self.min_seq_length = int(options[i])
			elif options[i] == '-isoutmsg': # Output messages to stdout(screen) or not. "0" means that PLEK be run quietly. "1" means that PLEK outputs the details of processing.
				i = i + 1
				self.isoutmsg = int(options[i])
			elif options[i] == '-isrmtempfile': # Remove temporary files or not. "0" means that PLEK retains temporary files. "1" means that PLEK remove temporary files.
				i = i + 1
				self.isrmtempfile = int(options[i])
			elif options[i] == '-nfold': # n-fold cross validation in search for optimal parameters.
				i = i + 1
				self.nfold = int(options[i]) 
			elif options[i] == '-k': # range of k. k=5 means that we will calculate usage of 1364 k-mer patterns. (k=1, 4 patterns; k=2, 16; k=3, 64; k=4, 256; k=5, 1024; 1364=4+64+256+1024)
				i = i + 1
				self.kmer = int(options[i])
			elif options[i] == '-model': # svm model file (output).
				i = i + 1
				self.modelfile = options[i]
			elif options[i] == '-range': # svm range file (output).
				i = i + 1
				self.rangefile = options[i]     
			elif options[i] == '-prefix': # The prefix of output files. 
				i = i + 1
				self.prefix = options[i]
			elif options[i] == '-isbalanced': # NOTE: isbalanced=1, need to be balanced. If the samples are unbalanced, it will subsample the overrepresented class to obtain an equal amount of positives and negatives;  isbalanced=0, not to be balanced. If the samples are unbalanced, svm will reweight the misclassification cost (adjust the wi parameters).
				i = i + 1
				self.is_mRNAlncRNA_balanced = int(options[i]) 
			elif options[i] == '-isrecompile': # re-compile source.
				i = i + 1
				self.is_recompile = int(options[i])
			elif options[i] == '-isadjustweight': # isadjustweight.
				i = i + 1
				self.isadjustweight = int(options[i])
			elif options[i] == '-x_fold': # 
				i = i + 1
				self.x_fold = int(options[i])				
			elif options[i] == '-islogging': # 
				i = i + 1
				self.islogging = int(options[i])	
			else:
				pass_through_opts.append(options[i])
			i = i + 1

		self.pass_through_string = ' '.join(pass_through_opts)


if __name__ == '__main__':
	def exit_with_help():
		print("""\
=====================
  USAGE AND EXAMPLES
=====================
Usage of PLEK Modelling -- used to build your own classifier with your 
                           training data (mRNA/lncRNA transcripts):
 
python PLEKModelling.py -mRNA mRNAs_fasta -lncRNA lncRNAs_fasta 
   -prefix prefix_of_output -log2c range_of_log2c -log2g range_of_log2g 
   -thread number_of_threads -model model_file -range range_file  
   -minlength min_length_of_sequence -isoutmsg 0_or_1 -isrmtempfile 0_or_1
   -k k_mer -nfold n_fold_cross_validation -isbalanced 0_or_1 
   
   -mRNA          mRNA transcripts used to build predictor, in fasta format.
   
   -lncRNA        lncRNA transcripts used to build predictor, in fasta format.
   
   -prefix        Prefix of the output files.
   
   -log2c        (Optional) The specified range of C parameter for the svm parameter 
                  search. Default value: 0,5,1. (from, to, by; 0,1,2,3,4,5)   
				  
   -log2g        (Optional) The specified range of G parameter for the svm parameter 
                  search. Default value: 0,-5,-1.(from, to, by; 0,-1,-2,-3,-4,-5) 
				  
   -thread       (Optional) The number of threads for running the PLEKModelling 
                  program. The bigger this number is, the faster PLEKModelling runs.
                  Note that a larger thread number means larger consumption of memory.
                  Default value: 12.
				  
   -model        (Optional) The name of a predictor model file (an output file
                  by PLEKModelling.py).   
				  
   -range        (Optional) The name of a svm range file (an output file by 
                  PLEKModelling.py).   
				  
   -minlength    (Optional) The minimum length of sequences. The sequences whose 
                  lengths are more than minlength will be processed. Default 
                  value: 200.             
				  
   -isoutmsg     (Optional) Output messages to stdout(screen) or not. "0" means 
                 that PLEKModelling be run quietly. "1" means that PLEKModelling 
                 outputs the details of processing. Default value: 0.   
				 
   -isrmtempfile (Optional) Remove temporary files or not. "0" means that PLEKModelling 
                  retains temporary files. "1" means that PLEKModelling remove temporary 
                  files. Default value: 1.
				  
   -k            (Optional) range of k. k=5 means that we will calculate usage of 
                 1364 k-mer patterns. (k=1, 4 patterns; k=2, 16; k=3, 64; k=4, 256; 
                 k=5, 1024; 1364=4+64+256+1024). Default value: 5. 
				 
   -nfold        (Optional) n-fold cross-validation in search for optimal parameters.
                 Default value: 10.   
   
   -isbalanced   (Optional) In the case of isbalanced=1, if the samples are 
                 unbalanced, it will subsample the overrepresented class to obtain an 
                 equal amount of positives and negatives.
                 Default value: 0.
           

Examples: 
1. $ python PLEKModelling.py -mRNA PLEK_mRNAs.fa -lncRNA PLEK_lncRNAs.fa -prefix 20140531 

   NOTE: To train a classifier using the mRNA sequences in the 'PLEK_mRNAs.fa' file
   and lncRNA in 'PLEK_lncRNAs.fa'. The program outputs the model in the file 
   '20140531.model' and the svm-scale range in '20140531.range'. 
   We can use the new model as follows:
    python PLEK.py -fasta PLEK_test.fa -out 20140531.predicted -thread 10  \\
    -range 20140531.range -model 20140531.model 
	
2. $ python PLEKModelling.py -mRNA PLEK_mRNAs.fa -lncRNA PLEK_lncRNAs.fa -prefix 20140601 \\
     -log2c 1,3,1 -log2g -1,-3,-1 -nfold 2 -k 4      

   NOTE: This example is used to demonstrate the usage of PLEKModelling.py
   in a simple/quick way. User can run this to check if our program can run correctly.
   To train a classifier using the mRNA sequences in the 'PLEK_mRNAs.fa' file
   and lncRNA in 'PLEK_lncRNAs.fa'. The range of log2c is 1,2,3. The range of log2g
   is -1,-2,-3. Use a 2-fold cross-validation. K is 4, it will calculate the usage
   of 340 patterns (4+16+63+256=340). The program outputs the model in the file 
   '20140601.model' and the svm-scale range in '20140601.range'. 
   We can use the new model as follows:
    python PLEK.py -fasta PLEK_test.fa -out 20140601.predicted -thread 10  \\
    -range 20140601.range -model 20140601.model -k 4

Notes:
   (1) In general, it is time-consuming to build a new classifier.          
   (2) The accuracy of a classifier model is connected with the quality and 
   quantity of training samples. We encourage users to supply as many reliable 
   samples as possible to PLEKModelling.py to train a classifier. We suggest
   that transcripts annotated with 'pseudogene', 'predicted' and 'putative' 
   be removed before training models. 
   (3) To parallel run PLEKModelling.py in PBS, please see 
   PLEK_howto_generate_scripts.pdf   
   
   
=====================
  CONTACTS
=====================
Aimin Li: LiAiminMail@gmail.com
Junying Zhang: jyzhang@mail.xidian.edu.cn     

======================
  WEBSITE
=====================
https://sourceforge.net/projects/plek/
           """)
		sys.exit(1)
		
   
	def compile_c(_opts): # re-compile source
		print('[{0}] Compiling svm, svm-train, svm-predict, svm-scale'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
		os.system("g++ -c " + _opts.script_dir + "svm.cpp -o " + _opts.script_dir + "svm.o")
		os.system("LNAG=C gcc -g -Wall " + _opts.script_dir + "svm-train.c " + _opts.script_dir + "svm.o -o " + _opts.script_dir + "svm-train -lstdc++ -lm")
		os.system("LNAG=C gcc -g -Wall " + _opts.script_dir + "svm-predict.c " + _opts.script_dir + "svm.o -o " + _opts.script_dir + "svm-predict  -lstdc++ -lm")
		os.system("LNAG=C gcc -g -Wall " + _opts.script_dir + "svm-scale.c " + _opts.script_dir + "svm.o -o " + _opts.script_dir + "svm-scale  -lstdc++ -lm")

		print('[{0}] Compiling PLEK_main, PLEK_spsn'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
		os.system("LNAG=C gcc -g -Wall " + _opts.script_dir + "PLEK_main.c -o " + _opts.script_dir + "PLEK -lm")
		os.system("LNAG=C gcc -g -Wall " + _opts.script_dir + "PLEK_spsn.c -o " + _opts.script_dir + "PLEK_spsn -lm ")		
	
	def sampling(pos_number, neg_number, prefix, x_fold):
		if pos_number> neg_number:
			os.system(' grep -e "^1"  '+prefix+'.scale | sort -R | sort -R | sed -n "1,'+str(x_fold*neg_number)+'p" > '+prefix+'.pos_selected ')
			os.system(' grep -e "^0"  '+prefix+'.scale  > '+prefix+'.neg_selected ')
			os.system(' cat '+prefix+'.pos_selected '+prefix+'.neg_selected > '+prefix+'.newscale ')
			os.system(' rm '+prefix+'.pos_selected ')
			os.system(' rm '+prefix+'.neg_selected ')
		if neg_number> pos_number:
			os.system(' grep -e "^0"  '+prefix+'.scale | sort -R | sort -R | sed -n "1,'+str(x_fold*pos_number)+'p" > '+prefix+'.neg_selected ')
			os.system(' grep -e "^1"  '+prefix+'.scale  > '+prefix+'.pos_selected ')
			os.system(' cat '+prefix+'.pos_selected '+prefix+'.neg_selected > '+prefix+'.newscale ')
			os.system(' rm '+prefix+'.pos_selected ')
			os.system(' rm '+prefix+'.neg_selected ')
			
	if len(sys.argv) < 2:
		exit_with_help()
	options = sys.argv
	try:

		print('[{0}] Beginning PLEK Modelling run (Version 1.2) '.format(time.strftime('%Y-%m-%d %H:%M:%S')))   
		
		# get input options
		_opts = GridOption(options);
		cmdline=None;
		
		c1,c2,step1=_opts.log2c.split(",")
		g1,g2,step2=_opts.log2g.split(",")
		grid_file_name="svm_grid_modelling.py"
		if c1==c2 and g1==g2:
			grid_file_name="svm_grid_modelling_singlet.py"
		
		logfile=_opts.prefix + "_modelling.logs"
		if _opts.islogging==1:
			if os.path.isfile(logfile):
				os.remove(logfile)
			if not os.path.isfile(logfile):
				open(logfile, 'a').close()
				with open(logfile, 'a') as log_file:
					log_file.write(time.strftime('%Y-%m-%d %H:%M:%S')+", started.\n")
		# recompile source
		if _opts.is_recompile or (not os.path.isfile(_opts.script_dir + 'PLEK')):
			compile_c(_opts)
		else:
			time.sleep(2)

		# check if the mRNA/lncRNA file exists
		if not os.path.isfile( _opts.mRNA_file):
			print("ERROR: No such file '" + _opts.mRNA_file + "'")
			with open(logfile, 'a') as log_file:
				log_file.write("ERROR: No such file '" + _opts.mRNA_file + "'\n")
			sys.exit(1)
		if not os.path.isfile( _opts.lncRNA_file):
			print("ERROR: No such file '" + _opts.lncRNA_file + "'")
			with open(logfile, 'a') as log_file:
				log_file.write("ERROR: No such file '" + _opts.lncRNA_file + "'\n")
			sys.exit(1) 
      
		if _opts.modelfile=="":
			_opts.modelfile=_opts.prefix + ".model"
		if _opts.rangefile=="":
			_opts.rangefile=_opts.prefix + ".range"  
	  
		# output options:
		print("	Settings:")
		print("	mRNAs:                        " + _opts.mRNA_file)
		print("	lncRNAs:                      " + _opts.lncRNA_file)
		print("	Minimal sequence length (nt): " + str(_opts.min_seq_length) )
		print("	Thread count:                 " + str(_opts.thread_count ))
		print("	n-fold cross validation:      " + str(_opts.nfold ))
		print("	log2c (from, to, by):         " + _opts.log2c )
		print("	log2g (from, to, by):         " + _opts.log2g )
		print("	k-mer:                        " + str(_opts.kmer ))
		print("	Prefix of output files:       " + _opts.prefix) 
		print("	Model (output):               " + _opts.modelfile )
		print("	Range (output):               " + _opts.rangefile )   
		with open(logfile, 'a') as log_file:
			log_file.write("	Settings:"+ "\n")
			log_file.write("	mRNAs:                        " + _opts.mRNA_file+ "\n")
			log_file.write("	lncRNAs:                      " + _opts.lncRNA_file+ "\n")
			log_file.write("	Minimal sequence length (nt): " + str(_opts.min_seq_length) + "\n")
			log_file.write("	Thread count:                 " + str(_opts.thread_count )+ "\n")
			log_file.write("	n-fold cross validation:      " + str(_opts.nfold )+ "\n")
			log_file.write("	log2c (from, to, by):         " + _opts.log2c + "\n")
			log_file.write("	log2g (from, to, by):         " + _opts.log2g + "\n")
			log_file.write("	k-mer:                        " + str(_opts.kmer )+ "\n")
			log_file.write("	Prefix of output files:       " + _opts.prefix+ "\n") 
			log_file.write("	Model (output):               " + _opts.modelfile + "\n")
			log_file.write("	Range (output):               " + _opts.rangefile + "\n") 
			
		svm_file=str(_opts.prefix)+"_allsvm";
		print('[{0}] PLEK Modelling 1.2.0 is running'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
		
		with open(logfile, 'a') as log_file:
			log_file.write(time.strftime('%Y-%m-%d %H:%M:%S')+", calculating k-mer usage frequencies.\n")
		
		# calculate k-mer usage frequencies, output 
		print('[{0}] Calculating k-mer usage'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
		if _opts.mRNA_file!="" and _opts.lncRNA_file!="" and _opts.is_mRNAlncRNA_balanced==1:  # balanced
			if not os.path.isfile(_opts.mRNA_file):
				print("ERROR: No such file '" + _opts.mRNA_file + "'")
				sys.exit(1)
			if not os.path.isfile(_opts.lncRNA_file):
				print("ERROR: No such file '" + _opts.lncRNA_file + "'")
				sys.exit(1)                
			cmdline = (_opts.script_dir + 'PLEK \\\n      -p {0} \\\n      -n {1} \\\n      -o {2} \\\n      -s 1 -d 5 -k {3} -l {4} -b -isoutmsg {5} -isrmtempfile {6}').format\
				(_opts.mRNA_file, _opts.lncRNA_file, _opts.prefix, _opts.kmer, _opts.min_seq_length, _opts.isoutmsg, _opts.isrmtempfile)
		if _opts.mRNA_file!="" and _opts.lncRNA_file!="" and _opts.is_mRNAlncRNA_balanced==0: # unbalanced
			if not os.path.isfile(_opts.mRNA_file):
				print("ERROR: No such file '" + _opts.mRNA_file + "'")
				sys.exit(1)
			if not os.path.isfile(_opts.lncRNA_file):
				print("ERROR: No such file '" + _opts.lncRNA_file + "'")
				sys.exit(1)              
			cmdline = (_opts.script_dir + 'PLEK \\\n      -p {0} \\\n      -n {1} \\\n      -o {2} \\\n      -s 1 -d 5 -k {3} -l {4}    -isoutmsg {5} -isrmtempfile {6}').format\
				(_opts.mRNA_file, _opts.lncRNA_file, _opts.prefix, _opts.kmer, _opts.min_seq_length, _opts.isoutmsg, _opts.isrmtempfile)

		if _opts.isoutmsg == 1:
			print("     " + cmdline)
		with open(logfile, 'a') as log_file:
			log_file.write("-------------\n\t"+cmdline+"\n-------------\n")
		os.system(cmdline) 


		
		# if unbalanced, get the number of positive and negative samples.
		# awk '$1==1' _allsvm | wc -l
		cmdline="awk '$1==1' " +  _opts.prefix + "_allsvm  | wc -l "
		pos_sample_no=int(commands.getstatusoutput(cmdline)[1])
		cmdline="awk '$1==0' " +  _opts.prefix + "_allsvm  | wc -l "
		neg_sample_no=int(commands.getstatusoutput(cmdline)[1])
		w1=1
		w0=1
		if pos_sample_no>300 and neg_sample_no>300: # cmd was executed successfully.
			if _opts.isadjustweight==1 and pos_sample_no>neg_sample_no:
				w0=pos_sample_no/neg_sample_no
			if _opts.isadjustweight==1 and pos_sample_no<neg_sample_no:
				w1=neg_sample_no/pos_sample_no
		else: # cmd failed, or all short sequences were removed , or 			
			print("Positive sample number: " + str(pos_sample_no)  )
			print("Negative sample number: " + str(neg_sample_no)  )
			print("WARNING: Number of samples is too small. " )
			sys.exit(1)

		with open(logfile, 'a') as log_file:
			log_file.write("\tPositive sample number: " +str(pos_sample_no) +";  Negative sample number: "  + str(neg_sample_no)+"\n")

		with open(logfile, 'a') as log_file:
			log_file.write(time.strftime('%Y-%m-%d %H:%M:%S')+", svm-scale.\n")
		
		# svm-scale
		print('[{0}] Scaling k-mer usage'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
		if not os.path.isfile( _opts.prefix + "_allsvm"):
			print("ERROR: No such file '" +  _opts.prefix + "_allsvm" + "'")
			sys.exit(1) 
		
		cmdline= (_opts.script_dir +'svm-scale \\\n      -l 0 -u 1 -s {0} \\\n      {1} \\\n      > {2} ').format\
				(_opts.rangefile, _opts.prefix + "_allsvm", _opts.prefix + ".scale")
		if _opts.isoutmsg == 1:
			print("     " + cmdline)
		with open(logfile, 'a') as log_file:
			log_file.write("-------------\n\t"+cmdline+"\n-------------\n")
		os.system(cmdline) 

		
		if not os.path.isfile( _opts.prefix + ".scale"):
			print("ERROR: No such file '" +  _opts.prefix + ".scale" + "'")
			sys.exit(1) 

		with open(logfile, 'a') as log_file:
			log_file.write(time.strftime('%Y-%m-%d %H:%M:%S')+", grid.py.\n")
		with open(logfile, 'a') as log_file:
			log_file.write("\tParameters is logged to " + _opts.prefix + ".scale.grid.out" + "\n")
			
		# grid.py, output to a file for subsequent analysis
		print('[{0}] Choosing optimal classifier parameters. This step is usally time-consuming'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
		
		new_pos_sample_no=0
		new_neg_sample_no=0
		 
		scaled_file=_opts.prefix + ".scale"
		if	pos_sample_no+neg_sample_no>50000 and (pos_sample_no>_opts.x_fold*neg_sample_no or pos_sample_no*_opts.x_fold<neg_sample_no ):
			sampling(pos_sample_no,neg_sample_no,_opts.prefix,_opts.x_fold)
			scaled_file=_opts.prefix + ".newscale"
			if pos_sample_no>neg_sample_no:
				new_pos_sample_no=_opts.x_fold*neg_sample_no
				new_neg_sample_no=neg_sample_no
			else:
				new_pos_sample_no=pos_sample_no
				new_neg_sample_no=_opts.x_fold*pos_sample_no
			
		memory=500+(new_pos_sample_no+new_neg_sample_no)*5000/100000
		epsilon =0.1
		if (new_pos_sample_no+new_neg_sample_no)<50000:
			epsilon =0.01
		if (new_pos_sample_no+new_neg_sample_no)<10000:
			epsilon =0.001			

		
		para_select_w1=1
		para_select_w0=1
		if _opts.isadjustweight==1 and new_pos_sample_no!=new_neg_sample_no:
			if new_pos_sample_no>new_neg_sample_no:
				para_select_w0=new_pos_sample_no/new_neg_sample_no
			else:
				para_select_w1=new_neg_sample_no/new_pos_sample_no

			
		if _opts.isadjustweight==1:
			cmdline= ("python " + _opts.script_dir +grid_file_name+" \\\n      -log2c {0} -log2g {1} -thread {2} \\\n      -svmtrain {3} \\\n      -gnuplot null -v {4} -e "+str(epsilon)+" -m "+str(memory)+" -l "+ logfile +" -w1 {5} -w0 {6} \\\n       {7}  ").format\
				(_opts.log2c, _opts.log2g, _opts.thread_count, _opts.svmtrain_pathname, _opts.nfold, para_select_w1, para_select_w0, scaled_file)
		else:
			cmdline= ("python " + _opts.script_dir +grid_file_name+" \\\n      -log2c {0} -log2g {1} -thread {2} \\\n      -svmtrain {3} \\\n      -gnuplot null -v {4} -e "+str(epsilon)+" -m "+str(memory)+" -l "+ logfile +" \\\n      {5}  ").format\
				(_opts.log2c, _opts.log2g, _opts.thread_count, _opts.svmtrain_pathname, _opts.nfold,  scaled_file)
				
		if _opts.isoutmsg == 1:
			print("     " + cmdline)
		with open(logfile, 'a') as log_file:
			log_file.write("-------------\n\t"+cmdline+"\n-------------\n")
		os.system(cmdline) 
		
		if os.path.isfile( _opts.prefix + ".newscale"):
			if os.path.isfile( _opts.prefix + ".scale.grid.out"):
				os.system(" rm " + _opts.prefix + ".scale.grid.out")
			if os.path.isfile( _opts.prefix + ".scale.grid.err"):
				os.system(" rm " + _opts.prefix + ".scale.grid.err")
			os.system(" mv " + _opts.prefix + ".newscale.grid.out " + _opts.prefix + ".scale.grid.out")
			os.system(" mv " + _opts.prefix + ".newscale.grid.err " + _opts.prefix + ".scale.grid.err")
       
		# get optimal parameter: c, g
		# cat 20140530.scale.grid.out | sed 's/parameters.*]//g' | sed 's/(best//g' | sed 's/c=//g' | sed 's/g=//g' | sed 's/rate=//g' | sed 's/)//g'| sed 's/,//g' | awk '{ print $3,$1,$2,$4,$5}'| sort -r | sed -n '1,1p' 
		# 91.4 1.0 -2.0 2.0 0.25
		cmdline="cat " + _opts.prefix + ".scale" +".grid.out | sed 's/parameters.*]//g' | sed 's/(best//g' | sed 's/c=//g' | sed 's/g=//g' | sed 's/rate=//g' | sed 's/)//g'| sed 's/,//g' | awk '{ print $3,$1,$2,$4,$5}' | sort -r | sed -n '1,1p' > " +  _opts.prefix + ".scale.grid.out.opt" 
		if _opts.isoutmsg == 1:
			print("     " + cmdline)
		with open(logfile, 'a') as log_file:
			log_file.write("-------------\n\t"+cmdline+"\n-------------\n")
		os.system(cmdline) 
		
		# output:  _opts.prefix + ".scale.grid.out.opt" 
		# rate, log2c, log2g, c, g
		rate=0.0
		log2c=0.0
		log2g=0.0
		c=0.0
		g=0.0
		for line in open(_opts.prefix + ".scale.grid.out.opt", 'r'):
			line = line.strip()  
			rate=float(line.split()[0])
			log2c=float(line.split()[1])
			log2g=float(line.split()[2])
			c=float(line.split()[3])
			g=float(line.split()[4])		
			

		with open(logfile, 'a') as log_file:
			log_file.write("\tOptimal parameters is logged to " + _opts.prefix + ".scale.grid.out.opt" + "\n")
			
		with open(logfile, 'a') as log_file:
			log_file.write(time.strftime('%Y-%m-%d %H:%M:%S')+", building model.\n")		
		
		if rate<=0:
			print "ERROR. rate<=0."
			sys.exit(1)
		
		print('[{0}] Building model'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
		
		# svm-train a model
		#./svm-train -c ${c} -g ${g} -m 3000  rout.scale xxx.model

		memory=500+(pos_sample_no+neg_sample_no)*5000/100000
		if _opts.isadjustweight==1: 
			cmdline= (_opts.script_dir +'svm-train  \\\n      -c {0} -g {1} -m '+str(memory)+' -w1 {2} -w0 {3} -q  \\\n      {4} \\\n      {5} ').format\
				   (c,g, w1, w0, _opts.prefix + ".scale", _opts.prefix + ".model")
		else:
			cmdline= (_opts.script_dir +'svm-train  \\\n      -c {0} -g {1} -m '+str(memory)+' -q  \\\n      {2} \\\n      {3} ').format\
				   (c,g, _opts.prefix + ".scale", _opts.prefix + ".model")
				
		if _opts.isoutmsg == 1:
			print("     " + cmdline)
		with open(logfile, 'a') as log_file:
			log_file.write("-------------\n\t"+cmdline+"\n-------------\n")
		if c<0 or g<0:
			print("\t"+"c<0 or g<0")
			sys.exit(1)
		os.system(cmdline)		
		
       
		# remove temporary files
		#  os.remove(_opts.prefix+'_temp_'+str(fn))
		if _opts.isrmtempfile==1:
			if os.path.isfile( _opts.prefix + "_allsvm"):		
				os.remove(_opts.prefix + "_allsvm")
			if os.path.isfile( _opts.prefix + "_allsvmdesc"):		
				os.remove(_opts.prefix + "_allsvmdesc")
			#if os.path.isfile( _opts.prefix + "_logs"):		
			#	os.remove(_opts.prefix + "_logs")
			if os.path.isfile( _opts.prefix + ".scale"):		
				os.remove(_opts.prefix + ".scale")
			if os.path.isfile( _opts.prefix + ".newscale"):		
				os.remove(_opts.prefix + ".newscale")				
			#if os.path.isfile( _opts.prefix + ".scale.grid.err"):		
			#	os.remove(_opts.prefix + ".scale.grid.err")
			#if os.path.isfile( _opts.prefix + ".scale.grid.out"):		
			#	os.remove(_opts.prefix + ".scale.grid.out")
			#if os.path.isfile( _opts.prefix + ".scale.grid.out.opt"):		
			#	os.remove(_opts.prefix + ".scale.grid.out.opt")
			if os.path.isfile( _opts.prefix + ".scale"):		
			    os.remove(_opts.prefix + ".scale")
			
		# output		 
		print('[{0}] Run complete'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
		print('	Model file: {0}'.format(_opts.modelfile))
		print('	Range file: {0}'.format(_opts.rangefile))
		
		with open(logfile, 'a') as log_file:
			log_file.write('	Model file: {0}'.format(_opts.modelfile) + "\n")
		with open(logfile, 'a') as log_file:
			log_file.write('	Range file: {0}'.format(_opts.rangefile) + "\n")

		with open(logfile, 'a') as log_file:
			log_file.write(time.strftime('%Y-%m-%d %H:%M:%S')+", end.\n")	
		
		sys.exit(0)
		
	except (IOError,ValueError) as e:
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write('Try "python {0}" for more information.\n'.format(sys.arg[0]))
		sys.exit(1)