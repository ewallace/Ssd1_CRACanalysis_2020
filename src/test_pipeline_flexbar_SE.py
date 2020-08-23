#!/usr/bin/python
__author__      = "Sander Granneman"
__copyright__   = "Copyright 2018"
__version__     = "0.5.1e"
__credits__     = ["Sander Granneman","Edward Wallace"]
__email__       = "sgrannem@staffmail.ed.ac.uk"
__status__      = "beta"

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

parser = cmdline.get_argparse(description="")
parser.add_argument("-f",dest="forwardreads",help="the path to your fastq read files.",metavar="data_1.fastq data_2.fastq ...",nargs="*",default=None)
parser.add_argument("--name",dest="name",help="provide a single word describing the run. Default is 'analysis' with a time stamp",default="analysis_%s" % time.strftime("%d%m%Y"))
parser.add_argument("-a","--adapter",dest="adapter",help="the path to the file containing the 3' adapter sequences for trimming the reads using flexbar. If you do not provide adapter sequences, the trimming step will be skipped",default=None)
parser.add_argument("-p","--processors",dest="processors",type=int,help="indicate how many processors you want to use for analyses. Default is 8",default=8)
args = parser.parse_args()


def runFlexBar(inputfile,outputfile):
    """ runs Flexbar on inputfile to remove the adapter sequence from the forward reads """
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

### start of pipeline

pipeline = Pipeline(name="Single-end data CRAC pipeline")
                

pipeline.transform(task_func = runFlexBar,
                    input  = startingfiles,
                    filter = formatter("^.+/([^/]+).(san)?fastq$"),                     
                    output = "%s/flexbar_trimmed/{1[0]}_trimmed.fastq" % root_dir,
                ).follows(pipeline.mkdir(os.path.join(root_dir,"flexbar_trimmed")))
                

pipeline_run(multiprocess=args.processors,verbose=5,checksum_level=3)
