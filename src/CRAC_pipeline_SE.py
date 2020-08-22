#!/usr/bin/python
__author__		= "Sander Granneman"
__copyright__	= "Copyright 2018"
__version__		= "0.5.0"
__credits__		= ["Sander Granneman"]
__email__		= "sgrannem@staffmail.ed.ac.uk"
__status__		= "beta"

from ruffus import *
from ruffus.cmdline import MESSAGE
from collections import defaultdict
import ruffus.cmdline as cmdline
import subprocess
import platform
import glob
import re
import os
import argparse
import sys
import time

DEFAULTADAPT = "TGGAATTCTCGGGTGCCAAGGC"

parser = cmdline.get_argparse(description="CRAC pipeline for processing single-end CRAC data")
parser.add_argument("-f",dest="forwardreads",help="the path to your fastq read files.",metavar="data_1.fastq data_2.fastq ...",nargs="*",default=None)
parser.add_argument("-g","--gtf",dest="gtf",help="the path to your gtf annotation file",metavar="rRNA.gtf",default=None)
parser.add_argument("-c","--chromosome",dest="chromosome",help="the path to your chromosome length file",metavar="chromosome_lengths.txt",default=None)
parser.add_argument("-i","--ignoremutations",dest="ignoremuts",action="store_true",help="to tell the novoalign parser to ignore mutations. Useful when analyzing ChemModSeq data. Default is OFF",default=False)
parser.add_argument("--nocollapse",dest="nocollaps",action="store_true",help="skips the trimming and collapsing step of the analysis. Note that files should have a .fasta file extension! Default is OFF",default=False)
parser.add_argument("--novoindex",dest="novoindex",help="the path to your novoindex file",metavar="yeast.novoindex",default=None)
parser.add_argument("--name",dest="name",help="provide a single word describing the run. Default is 'analysis' with a time stamp",default="analysis_%s" % time.strftime("%d%m%Y"))
parser.add_argument("--sgr",dest="sgr",help="to make sgr files with read counts for each genomic position. Default is off",action="store_true",default=False)
parser.add_argument("-b","--barcodes",dest="barcodes",help="the path to the file containing the list of barcodes. If you do not provide a barcode file, the demultiplexing step will be skipped",metavar="barcodes.txt",default=None)
parser.add_argument("-a","--adapter",dest="adapter",help="the path to the file containing the 3' adapter sequences for trimming the reads using flexbar. If you do not provide adapter sequences, the trimming step will be skipped",default=None)
parser.add_argument("--truseq",dest="truseq",action="store_true",help="add this flag if your library was prepared using TruSeq kits. NOTE! This requires Flexbar version 3.4.0 or later! Default is False",default=False)
parser.add_argument("-m","--mismatches",dest="mismatches",type=int,help="indicate how many mismatches you allow for demultiplexing. Default is 1",default=1)
parser.add_argument("-p","--processors",dest="processors",type=int,help="indicate how many processors you want to use for analyses. Default is 8",default=8)
args = parser.parse_args()

def getFileBaseName(filename):
	""" removes path and extension from file name """
	return os.path.splitext(os.path.basename(filename))[0]

def getBarcodeInfo(barcodefile):
	return ["%s.fastq" % "_".join(line.strip().split()) for line in open(barcodefile,"r").readlines()]
	
def runFlexBar(inputfile,outputfile):
	""" runs Flexbar on the data to remove the adapter sequence from the forward reads """
	outputfilename = "%s/%s_trimmed" % (os.path.join(root_dir,"flexbar_trimmed"),re.search("^.+/([^/]+).(san)?fastq(.gz)?",inputfile).group(1))
	if args.adapter:
		cmd = "flexbar -r '%s' -qf i1.8 -t '%s' -n 10 -ao 7 -a '%s' -qt 30" % (inputfile,outputfilename,args.adapter)
	elif args.truseq:
		cmd = "flexbar -r '%s' -qf i1.8 -t '%s' -n 10 -ao 7 -aa TruSeq -qt 30" % (inputfile,outputfilename)
		print cmd
	else:
		cmd = "flexbar -r '%s' -qf i1.8 -t '%s' -n 10 -ao 7 -a '%s' -qt 30" % (inputfile, outputfilename)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s" % outputfile)

def demultiplexSamples(inputfile,outputfile):
	""" demultiplexes all the samples """
	os.chdir(os.path.join(root_dir,"demultiplexed"))
	cmd = "pyBarcodeFilter.py -f '%s' -b '%s' -m '%s'" % (inputfile,os.path.join(home_dir,args.barcodes),args.mismatches)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s" % outputfile)
	os.chdir(root_dir)

def trimAndCollapse(inputfile,outputfile):
	""" Removes the last random nucleotide from the forward read, collapses the data and then splits it again into fasta files """
	cmd = "pyFastqDuplicateRemover.py -f '%s' -o '%s'" % (inputfile,outputfile)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s" % outputfile)

	
def alignReads(inputfile,outputfile):
	""" Runs novoalign on all the collapsed files"""
	cmd = "novoalign -d '%s' -f '%s' -r Random > '%s'" % (args.novoindex,inputfile,outputfile)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s" % outputfile)

def runPyReadCounters(inputfile,outputfile):
	""" runs pyReadCounters on all the files """
	outputfile = re.search("(.*)_count_output_reads.gtf",outputfile).group(1)
	if args.ignoremuts:
		string = "--mutations=nomuts"
	else:
		string = ""
	cmd = "pyReadCounters.py -f '%s' --gtf '%s' -v --rpkm -o '%s' '%s'" % (inputfile,args.gtf,outputfile,string)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s_count_output_reads.gtf" % outputfile)
	
def runSecondPyReadCounters(inputfile,outputfile):
	""" runs pyReadCounters on all the files """
	outputfile = re.search("(.*)_count_output_cDNAs.gtf",outputfile).group(1)
	cmd = "pyReadCounters.py -f '%s' --gtf '%s' -v --rpkm -o '%s' --mutations  nomuts --blocks" % (inputfile,args.gtf,outputfile)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s_count_output_reads.gtf" % outputfile)
	
def makeCoverageSgrFiles(inputfile,outputfiles):
	""" takes the pyReadCounters gtf output files and generates coverage sgr files """
	outputfile = re.search("(.*)_plus_strand_reads.sgr",outputfiles[0]).group(1)
	cmd = "pyGTF2sgr.py --gtf '%s' --zeros --count -v -o '%s' -c '%s'" % (inputfile,outputfile,args.chromosome)
	logger.info(cmd)
	os.system(cmd)
	#for i in outputfiles:
	#	os.system("touch %s" % i)
		
def makeCoverageBedgraphFiles(inputfile,outputfiles):
	""" takes the pyReadCounters gtf output files and generates bedgraph files for viewing of the data in genome browsers """
	outputfile = re.search("(.*)_plus_strand_reads.bedgraph",outputfiles[0]).group(1)
	cmd = "pyGTF2bedGraph.py --gtf '%s' --count -v --permillion -o '%s' -c '%s'" % (inputfile,outputfile,args.chromosome)
	logger.info(cmd)
	os.system(cmd)
	#for i in outputfiles:
	#	os.system("touch %s" % i)

### checking if dependencies are installed

dependencies = ["flexbar","pyGTF2sgr.py","pyGTF2bedGraph.py","pyReadCounters.py","pyBarcodeFilter.py","pyFastqJoiner.py","pyFastqSplitter.py","TrimNucs.py","novoalign","pyFastqDuplicateRemover.py"]

for i in dependencies:
	cmd = "where" if platform.system() == "Windows" else "which"
	try: 
		subprocess.call([cmd, "%s"])
	except: 
		sys.stderr.write("%s is not installed (properly) on your system. Please (re)install" % i)
		exit()
		
### start of pipeline

# setting the main working directory
home_dir = os.getcwd()
root_dir = "%s/%s" % (home_dir,args.name)
try:
	os.mkdir(root_dir)
except OSError:
	pass

### setting up logging
if not args.log_file:
	logfile = "%s/%s" % (root_dir,"log_file.txt")
else:
	logfile = "%s/%s" % (root_dir,args.log_file)
logger, logger_mutex = cmdline.setup_logging ("ChemModSeqPipeline",logfile,10)

### setting the starting files

startingfiles = [os.path.abspath(i) for i in args.forwardreads]

### checking if the correct pyBarcodeFilter version is installed. Need version 3.0 or higher!

result = subprocess.Popen(['pyBarcodeFilter.py', '--version'], stdout=subprocess.PIPE)
out,err = result.communicate()
print float(out)
if float(out) < 3.0:
	sys.stderr.write("To run this script you need to have pyBarcodeFilter version 3.0 or later installed. Please install the latest version of pyCRAC\n")
	exit()	

### start of pipeline

pipeline = Pipeline(name="Single-end data CRAC pipeline")
				

pipeline.transform(task_func = runFlexBar,
					input  = startingfiles,
					filter = formatter("^.+/([^/]+).(san)?fastq$"),						
					output = "%s/flexbar_trimmed/{1[0]}_trimmed.fastq" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"flexbar_trimmed")))
				
startingfiles = runFlexBar

if args.barcodes:
	pipeline.subdivide(task_func = demultiplexSamples,
						input  = startingfiles,
						filter = formatter("^.+/([^/]+).fastq"),
						output = ["%s/demultiplexed/{1[0]}_%s" % (root_dir,barcodestring) for barcodestring in getBarcodeInfo(args.barcodes)]
					).follows(pipeline.mkdir(os.path.join(root_dir,"demultiplexed")))
	startingfiles = demultiplexSamples
	
if not args.nocollaps:
	pipeline.transform(task_func = trimAndCollapse,
						input  = startingfiles,
						filter = formatter("^.+/([^/]+)_trimmed(?P<BARCODE>.*).fastq"),						
						output = "%s/collapsed/{1[0]}_trimmed{BARCODE[0]}.fasta" % root_dir
					).follows(pipeline.mkdir(os.path.join(root_dir,"collapsed")))
	
	startingfiles = trimAndCollapse

pipeline.transform( task_func = alignReads,
					input  = startingfiles,
					filter = formatter("^.+/([^/]+).fasta"),
					output = "%s/novo_files/{1[0]}.novo" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"novo_files")))

pipeline.transform( task_func = runPyReadCounters,
					input  = alignReads,
					filter = formatter("^.+/([^/]+).novo"),
					output = "%s/pyReadCounters_analyses/{1[0]}_count_output_reads.gtf" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"pyReadCounters_analyses")))

pipeline.transform(task_func = runSecondPyReadCounters,
					input  = alignReads,
					filter = formatter("^.+/([^/]+).novo"),
					output = "%s/pyReadCounters_blocks_nomuts_analyses/{1[0]}_blocks_nomuts_count_output_cDNAs.gtf" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"pyReadCounters_blocks_nomuts_analyses"))).follows(alignReads)

pipeline.subdivide(task_func = makeCoverageBedgraphFiles,
					input  = runSecondPyReadCounters,
					filter = formatter("^.+/([^/]+)_count_output_cDNAs.gtf"),
					output = ["%s/bedgraph_files/{1[0]}_plus_strand_reads.bedgraph" % root_dir,"%s/bedgraph_files/{1[0]}_minus_strand_reads.bedgraph" % root_dir]
				).follows(pipeline.mkdir(os.path.join(root_dir,"bedgraph_files"))).follows(runSecondPyReadCounters)

if args.sgr:
	pipeline.subdivide(task_func = makeCoverageSgrFiles,
						input  = runSecondPyReadCounters,
						filter = formatter("^.+/([^/]+)_count_output_cDNAs.gtf"),
						output = ["%s/sgr_files/{1[0]}_plus_strand_reads.sgr" % root_dir,"%s/sgr_files/{1[0]}_minus_strand_reads.sgr" % root_dir]
					).follows(pipeline.mkdir(os.path.join(root_dir,"sgr_files"))).follows(runSecondPyReadCounters)				

pipeline_run(multiprocess=args.processors,verbose=5,checksum_level=3)
