#!/usr/bin/env python

from optparse import OptionParser
import os

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(False)

import sys
sys.path.append("../../utils/")

# vmmana
import cluster as vmmcluster
import hit as vmmhit
import utils

ana_name = "t6t7ana"
def info(output) :
    print "%s    %s"%(ana_name, output)

if __name__=="__main__" :
    print "starting..."

    parser = OptionParser()
    parser.add_option("-i", "--input", default="")
    parser.add_option("-d", "--debug", action="store_true", default=False)
    parser.add_option("-n", "--nEvents", default=-1)
    parser.add_option("-c", "--calibFile", default="")
    (options, args) = parser.parse_args()

    input_file = options.input 
    dbg = options.debug
    nevents = options.nEvents
    calibFile = options.calibFile

    if not utils.file_ok(input_file) :
        info("exiting...")
        sys.exit()

