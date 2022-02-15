#!/usr/bin/python
__author__		= "Sander Granneman"
__copyright__	= "Copyright 2020"
__version__		= "0.7.0"
__credits__		= ["Sander Granneman", "Edward Wallace"]
__email__		= ["sgrannem@staffmail.ed.ac.uk", "edward.wallace@ed.ac.uk"]
__status__		= "beta"

""" CRAC_pipeline_SE_demult_dedup.py

This is a pipeline for processing single-end multiplexed CRAC sequencing data.
It is forked from CRAC_pipeline_SE.py v0.6.0, also written in ruffus.
This pipeline is specialised for single-end data that is multiplexed in-line, so
requires trimming of 3' adapters and demultiplexing, including deduplication. This 
simplifies the logic of the pipeline, no conditional execution needed.

"""

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

GFF_FIELDS = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]

parser = cmdline.get_argparse(description="CRAC pipeline for processing single-end multiplexed CRAC data")
parser.add_argument("-f",dest="forwardreads",help="the path to your fastq read files.",metavar="data_1.fastq data_2.fastq ...",nargs="*",default=None)
parser.add_argument("-g","--gtf",dest="gtf",help="the path to your gtf annotation file",metavar="rRNA.gtf",default=None)
parser.add_argument("-c","--chromosome",dest="chromosome",help="the path to your chromosome length file",metavar="chromosome_lengths.txt",default=None)
parser.add_argument("--novoindex",dest="novoindex",help="the path to your novoindex file",metavar="yeast.novoindex",default=None)
parser.add_argument("--name",dest="name",help="provide a single word describing the run. Default is 'analysis' with a time stamp",default="analysis_%s" % time.strftime("%d%m%Y"))
parser.add_argument("--sgr",dest="sgr",help="to make sgr files with read counts for each genomic position. Default is off",action="store_true",default=False)
parser.add_argument("-b","--barcodes",dest="barcodes",help="the path to the file containing the list of barcodes. If you do not provide a barcode file, the demultiplexing step will fail",metavar="barcodes.txt",default=None)
parser.add_argument("-a","--adapterfile",dest="adapterfile",help="file containing the 3' adapter sequences for trimming the reads using flexbar. If you do not provide adapter file or preset",default=None)
parser.add_argument("-aa","--adapterpreset",dest="adapterpreset",help="adapter preset string [TruSeq, SmallRNA, Methyl, Ribo, Nextera, and NexteraMP]. Requires Flexbar version 3.4.0 or later.",default=None)
parser.add_argument("-m","--mismatches",dest="mismatches",type=int,help="indicate how many mismatches you allow for demultiplexing. Default is 1",default=1)
parser.add_argument("-p","--processors",dest="processors",type=int,help="indicate how many processors you want to use for analyses. Default is 8",default=8)
parser.add_argument("--genelist",dest="genelist",help="provide a genelist file for pileup analysis, with one gene name per line,",default=None)
parser.add_argument("--genometab",dest="genometab",help="provide a genome sequence file in pyCRAC tab format for pileup analysis",default=None)
parser.add_argument("--transcriptgff",dest="transcriptgff",help="provide a gff file containing only transcript sequences for count analysis",default=None)
args = parser.parse_args()

def getBarcodeInfo(barcodefile):
	return ["%s.fastq" % "_".join(line.strip().split()) for line in open(barcodefile,"r").readlines()]
	
def runFlexBar(inputfile,outputfile,qflags = "-q TAIL -qf i1.8 -qt 30 --qtrim-post-removal"):
	""" runs Flexbar on the data to remove the adapter sequence from the forward reads """
	if args.adapterfile:
		cmd = "flexbar -r '%s' --output-reads '%s' -n 10 %s -ao 5 --adapters '%s'" % (inputfile,outputfile,qflags,args.adapterfile)
	elif args.adapterpreset:
		cmd = "flexbar -r '%s' --output-reads '%s' -n 10 '%s' -ao 5 -aa '%s'" % (inputfile,outputfile,qflags,args.adapterpreset)
	else:
		cmd = "flexbar -r '%s' --output-reads '%s' -n 10 '%s' -ao 5" % (inputfile,outputfile,qflags)
	logger.info(cmd)
	os.system(cmd)

def demultiplexSamples(inputfile,outputfile):
	""" demultiplexes all the samples """
	os.chdir(os.path.join(root_dir,"demultiplexed"))
	cmd = "pyBarcodeFilter.py -f '%s' -b '%s' -m '%s'" % (inputfile,os.path.join(home_dir,args.barcodes),args.mismatches)
	logger.info(cmd)
	os.system(cmd)
	os.chdir(root_dir)

def collapseDuplicates(inputfile,outputfile):
	""" Removes the last random nucleotide from the forward read, collapses the data and then splits it again into fasta files """
	cmd = "pyFastqDuplicateRemover.py -f '%s' -o '%s'" % (inputfile,outputfile)
	logger.info(cmd)
	os.system(cmd)
	
def alignReads(inputfile,outputfile,alignfiletype = "SAM"):
	""" Runs novoalign on all the collapsed files"""
	cmd = "novoalign -c 1 -d '%s' -o '%s' -f '%s' -r Random > '%s'" % (args.novoindex,alignfiletype,inputfile,outputfile)
	logger.info(cmd)
	os.system(cmd)

def runPyReadCounters(inputfile,outputfile,alignfiletype="sam"):
	""" runs pyReadCounters on all the files """
	outputfile = re.search("(.*)_count_output_reads.gtf",outputfile).group(1)
	cmd = "pyReadCounters.py -f '%s' --gtf '%s' -v --rpkm -o '%s' --file_type='%s' " % (inputfile,args.gtf,outputfile,alignfiletype)
	logger.info(cmd)
	os.system(cmd)
	
def runPyReadCountersBlocksNoMuts(inputfile,outputfile,alignfiletype="sam"):
	""" runs pyReadCounters on all the files with options --mutations  nomuts --blocks """
	outputfile = re.search("(.*)_count_output_cDNAs.gtf",outputfile).group(1)
	cmd = "pyReadCounters.py -f '%s' --gtf '%s' -v --rpkm -o '%s' --file_type='%s' --mutations  nomuts --blocks" % (inputfile,args.gtf,outputfile,alignfiletype)
	logger.info(cmd)
	os.system(cmd)
	
def makeCoverageSgrFiles(inputfile,outputfiles):
	""" takes the pyReadCounters gtf output files and generates coverage sgr files """
	outputfile = re.search("(.*)_plus.sgr",outputfiles[0]).group(1)
	cmd = "pyGTF2sgr.py --gtf '%s' --zeros --count -v -o '%s' -c '%s'" % (inputfile,outputfile,args.chromosome)
	logger.info(cmd)
	os.system(cmd)
		
def makeCoverageBedgraphFiles(inputfile,outputfiles):
	""" takes the pyReadCounters gtf output files and generates bedgraph files for viewing of the data in genome browsers """
	outputfile = re.search("(.*)_plus.bedgraph",outputfiles[0]).group(1)
	cmd = "pyGTF2bedGraph.py --gtf '%s' --count -v --permillion -o '%s' -c '%s'" % (inputfile,outputfile,args.chromosome)
	logger.info(cmd)
	os.system(cmd)

def runPyCalculateFDRs(inputfile,outputfile):
	""" runs pyCalculateFDRs.py to find peaks with false discovery rates, on protein-coding genes """
	cmd = "pyCalculateFDRs.py -f '%s' -o '%s' -c '%s' --gtf '%s' -a protein_coding -m 0.01 --min=5 " % (inputfile,outputfile,args.chromosome, args.gtf)
	logger.info(cmd)
	os.system(cmd)

def runPyPileup(inputfile,outputfile,alignfiletype = "sam"):
	""" runs pyPileup.py to make pileup tables of reads and deletions/mutations for a given genelist """
	cmd = "pyPileup.py -f '%s' -o '%s' --tab '%s' --gtf '%s' --file_type='%s' -g '%s'" % (inputfile,outputfile,args.genometab, args.gtf, alignfiletype, args.genelist)
	logger.info(cmd)
	os.system(cmd)
	
def runHTSeqCountTranscripts(inputfile,outputfile,minaqual=2):
	""" runs htseq-count to assign reads to primary transcript features """
	# Note: this did not play well with novoalign reads, so not being run for now.
	cmd = "htseq-count -n 8 -a '%s' --counts_output '%s' '%s' '%s'  " % (minaqual, outputfile, inputfile, args.gtf)
	logger.info(cmd)
	os.system(cmd)

def sortBamFiles(inputfile,outputfile):
	""" sorts .sam file, to produce a sorted compressed .bam """
	outputfile = outputfile
	cmdsort = "samtools view -b '%s' | samtools sort -@ 3 -O bam -o '%s' -T '%s'/tmp - " % (inputfile, outputfile, root_dir)
	logger.info(cmdsort)
	os.system(cmdsort)

def indexBamFiles(inputfile,outputfile):
	""" indexes a sorted .bam file, to produce a matching .bam.bai """
	cmdindex = "samtools index '%s' " % (inputfile)
	logger.info(cmdindex)
	os.system(cmdindex)

def runFastQC(inputfile,outputfile):
	""" runs fastQC on input fastq file """
	outputdir = os.path.split(outputfile)[0]
	cmd = "fastqc -o '%s' '%s'" % (outputdir, inputfile)
	logger.info(cmd)
	os.system(cmd)

def runBamQC(inputfile,outputfile):
	""" runs BamQC on input bam/sam file """
	outputdir = os.path.split(outputfile)[0]
	cmd = "bamqc -o '%s' '%s'" % (outputdir, inputfile)
	logger.info(cmd)
	os.system(cmd)

def runMultiCovTranscript(inputfiles,outputfile):
	""" counts reads to transcript features with bedtools multiBamCov """
	inputfileorder = ' '.join(inputfiles)
	inputfilebasenames = [ os.path.basename(f) for f in inputfiles ] 
	logger.info("multiBamCov file order is %s" % ' '.join(inputfilebasenames))
	with open(outputfile, "w") as outhandle:
		outhandle.write( "# allsample_transcriptcounts.txt - counts to transcripts from multiBamCov\n")
		outhandle.write( "# Run from CRAC_pipeline_SE_demult_dedup.py\n")
		outhandle.write( "# %s, date %s\n" % (args.name, time.strftime("%d-%m-%Y")))
		outhandle.write( "%s\n" % '\t'.join(GFF_FIELDS + inputfilebasenames))
	cmd = "multiBamCov -s -bed '%s' -bams %s >> '%s'" % (args.transcriptgff, inputfileorder, outputfile)
	logger.info(cmd)
	os.system(cmd)

def makeGenomeCoverageBedgraph(inputfile,outputfiles):
	""" takes the sorted indexed bam files and generates bedgraph files for viewing of the data in genome browsers """
	outputfile = re.search("(.*)_plus.bedgraph",outputfiles[0]).group(1)
	cmdplus = "genomeCoverageBed -bg -strand + -ibam %s > %s" % (inputfile,outputfiles[0])
	logger.info(cmdplus)
	os.system(cmdplus)
	cmdminus = "genomeCoverageBed -bg -strand - -ibam %s > %s" % (inputfile,outputfiles[1])
	logger.info(cmdminus)
	os.system(cmdminus)

### checking if dependencies are installed

dependencies = ["flexbar","pyGTF2sgr.py","pyGTF2bedGraph.py","pyReadCounters.py","pyBarcodeFilter.py","pyFastqJoiner.py","pyFastqSplitter.py","novoalign","pyFastqDuplicateRemover.py","samtools","fastqc","bamqc","multiBamCov"]

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
logger, logger_mutex = cmdline.setup_logging ("CRAC_Pipeline_SE",logfile,10)

### setting the starting files

startingfiles = [os.path.abspath(i) for i in args.forwardreads]

logger.info("analysing the following files:\n%s\n" % "\n".join(startingfiles))

### start of pipeline definition

pipeline = Pipeline(name="Single-end data CRAC pipeline for multiplexed libraries")

logger.info("Trimming the reads")
pipeline.transform(task_func = runFlexBar,
					input  = startingfiles,
					filter = formatter("^.+/([^/]+).(san)?fastq"),						
					output = "%s/flexbar_trimmed/{1[0]}_trimmed.fastq" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"flexbar_trimmed")))
			
logger.info("Demultiplexing the reads")
pipeline.subdivide(task_func = demultiplexSamples,
					input  = runFlexBar,
					filter = formatter("^.+/([^/]+).fastq"),
					output = ["%s/demultiplexed/{1[0]}_%s" % (root_dir,barcodestring) for barcodestring in getBarcodeInfo(args.barcodes)]
				).follows(pipeline.mkdir(os.path.join(root_dir,"demultiplexed")))

#logger.info("demultiplexed fastq quality control")
pipeline.transform( task_func = runFastQC,
					input  = demultiplexSamples,
					filter = formatter("^.+/([^/]+).fastq"),
					output = "%s/demultiplexed_fastqc/{1[0]}_fastqc.html" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"demultiplexed_fastqc")))

pipeline.transform(task_func = collapseDuplicates,
					input  = demultiplexSamples,
					filter = formatter("^.+/([^/]+)_trimmed(?P<BARCODE>.*).fastq"),						
					output = "%s/collapsed/{1[0]}_trimmed{BARCODE[0]}.fasta" % root_dir
				).follows(pipeline.mkdir(os.path.join(root_dir,"collapsed")))

logger.info("Aligning reads")
pipeline.transform( task_func = alignReads,
					input  = collapseDuplicates,
					filter = formatter("^.+/([^/]+).fasta"),
					output = "%s/aligned_sam/{1[0]}.sam" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"aligned_sam")))

#logger.info("aligned reads quality control")
pipeline.transform( task_func = runBamQC,
					input  = alignReads,
					filter = formatter("^.+/([^/]+).sam"),
					output = "%s/aligned_bamqc/{1[0]}_bamqc.html" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"aligned_bamqc")))

pipeline.transform( task_func = sortBamFiles,
					input  = alignReads,
					filter = formatter("^.+/([^/]+).sam"),
					output = "%s/aligned_bamsorted/{1[0]}.bam" % root_dir, 
				).follows(pipeline.mkdir(os.path.join(root_dir,"aligned_bamsorted")))

pipeline.transform( task_func = indexBamFiles,
					input  = sortBamFiles,
					filter = formatter("^.+/([^/]+).bam"),
					output = "%s/aligned_bamsorted/{1[0]}.bam.bai" % root_dir, 
				)

logger.info("Mapping reads to genomic features")
pipeline.transform( task_func = runPyReadCounters,
					input  = alignReads,
					filter = formatter("^.+/([^/]+).sam"),
					output = "%s/pyReadCounters_analyses/{1[0]}_count_output_reads.gtf" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"pyReadCounters_analyses"))).follows(alignReads)

pipeline.transform(task_func = runPyReadCountersBlocksNoMuts,
					input  = alignReads,
					filter = formatter("^.+/([^/]+).sam"),
					output = "%s/pyReadCountersBlocksNoMuts_analyses/{1[0]}_blocks_nomuts_count_output_cDNAs.gtf" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"pyReadCountersBlocksNoMuts_analyses"))).follows(alignReads)

pipeline.subdivide(task_func = makeCoverageBedgraphFiles,
					input  = runPyReadCountersBlocksNoMuts,
					filter = formatter("^.+/([^/]+)_count_output_cDNAs.gtf"),
					output = ["%s/bedgraph_files/{1[0]}_plus.bedgraph" % root_dir,"%s/bedgraph_files/{1[0]}_minus.bedgraph" % root_dir]
				).follows(pipeline.mkdir(os.path.join(root_dir,"bedgraph_files"))).follows(runPyReadCountersBlocksNoMuts)

logger.info("Making bedgraph files")

## Add: pyPileUp.py to calculate pileups with and without mutations/deletions on select genes
if args.genelist:
	pipeline.transform(task_func = runPyPileup,
						input  = alignReads,
						filter = formatter("^.+/([^/]+).sam"),
						output = "%s/pyPileup_analyses/{1[0]}_pileups.txt" % root_dir,
					).follows(pipeline.mkdir(os.path.join(root_dir,"pyPileup_analyses"))).follows(alignReads)

## Add: pyCalculateFDRs.py to calculate peaks and their occupancy
pipeline.subdivide(task_func = runPyCalculateFDRs,
					input  = runPyReadCounters,
					filter = formatter("^.+/([^/]+)_count_output_reads.gtf"),
					output = "%s/pyCalculateFDRs_analyses/{1[0]}_output_FDRs.gtf" % root_dir
				).follows(pipeline.mkdir(os.path.join(root_dir,"pyCalculateFDRs_analyses"))).follows(runPyReadCounters)

## Add: something to add up the number of reads per feature. Maybe genes inc. UTRs. 
# htseq-count was giving weird results with most reads unaligned, and aligned reads not mapped to features.
# placed on branch EW_SE_htseq and abandoned for now
# pipeline.transform(task_func = runHTSeqCountTranscripts,
# 					input  = alignReads,
# 					filter = formatter("^.+/([^/]+).sam"),
# 					output = "%s/HTSeqCountTranscript_analysis/{1[0]}_transcript_counts.txt" % root_dir,
# 				).follows(pipeline.mkdir(os.path.join(root_dir,"HTSeqCountTranscript_analysis"))).follows(alignReads)

if args.transcriptgff:
	pipeline.merge(task_func = runMultiCovTranscript,
						input  = sortBamFiles,
						output = "%s/multicov_analyses/allsample_transcriptcounts.txt" % root_dir,
					).follows(pipeline.mkdir(os.path.join(root_dir,"multicov_analyses"))).follows(indexBamFiles)

pipeline.subdivide(task_func = makeGenomeCoverageBedgraph,
					input  = sortBamFiles,
					filter = formatter("^.+/([^/]+).bam"),
					output = ["%s/bedgraph_genomecov/{1[0]}_plus.bedgraph" % root_dir,"%s/bedgraph_genomecov/{1[0]}_minus.bedgraph" % root_dir]
				).follows(pipeline.mkdir(os.path.join(root_dir,"bedgraph_genomecov"))).follows(indexBamFiles)

### print out flowchart describing pipeline
pipeline_printout_graph("%s/flowchart.png" %root_dir)

### run the pipeline
pipeline_run(multiprocess=args.processors,verbose=5,checksum_level=3)
