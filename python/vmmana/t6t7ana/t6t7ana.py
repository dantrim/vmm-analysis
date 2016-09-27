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

###############################################################################
# container class for holding event information
class AnalysisObject :
    def __init__(self, name_) :
        self.name = name_
        self.dbg = False

        self.run_number = -1

        self.do_tdo_calib = False
        self.tdo_file = ""
        self.do_pdo_calib = False
        self.pdo_file = ""

        self.tree_file = None
        self.tree = None

        # container for all hits
        self.hits = []
    
        # container for all clusters
        self.clusters = []

        # event data
        self.chip = []

    def load_event(self, event) :
        for ch in event.chip :
            self.chip.append(ch)

    def clear_containers(self) :
        self.hits = []
        self.clusters = []

    def load_tree(self, treename="", filename="") :
        rfile = r.TFile(filename)
        self.tree_file = rfile
        chain = r.TChain(treename)
        chain.Add(filename)
        self.tree = chain

    def load_tdo_calibration(self, do_calib=False, calib_file="") :
        print "load_tdo_calibration    NOT YET IMPLEMENTED"
        sys.exit()

    def load_pdo_calibration(self, do_calib=False, calib_file="") :
        print "load_pdo_calibration    NOT YET IMPLEMENTED"
        sys.exit()
        
###############################################################################
# event loop
def process_event(vmm, event_no, vmmdata) :
    info("process_event    %d"%event_no)
    vmm.load_event(vmmdata)
    print vmm.chip
    sys.exit()
    
    


###############################################################################
if __name__=="__main__" :
    print "starting..."

    parser = OptionParser()
    parser.add_option("-i", "--input", default="")
    parser.add_option("-d", "--debug", action="store_true", default=False)
    parser.add_option("-n", "--nEvents", default=-1)
    parser.add_option("--tdoCalibFile", default="")
    parser.add_option("--pdoCalibFile", default="")
    (options, args) = parser.parse_args()

    input_file = options.input 
    dbg = options.debug
    nevents = int(options.nEvents)
    tdoCalibFile = options.tdoCalibFile
    pdoCalibFile = options.pdoCalibFile

    if not utils.file_ok(input_file) :
        info("exiting...")
        sys.exit()

    # global analysis object
    vmmana = AnalysisObject("t6t7ana")
    vmmana.dbg = dbg

    do_tdo_calibration = False
    if tdoCalibFile!="" and not utils.file_ok(tdoCalibFile) :
        info("Will not apply tdo calibration constants")
    elif tdoCalibFile!="" and utils.file_ok(tdoCalibFile) :
        info("Will apply tod calibration")
        do_tdo_calibration = True

    do_pdo_calibration = False
    if pdoCalibFile!="" and not utils.file_ok(pdoCalibFile) :
        info("Will not apply pdo calibration constants")
    elif pdoCalibFile!="" and utils.file_ok(pdoCalibFile) :
        info("Will apply pdo calibration")
        do_pdo_calibration = True

    info("TODO: implement loading of calibration constants")
    #vmmana.load_tdo_calibration(do_tdo_calibration, tdoCalibFile)
    #vmmana.load_pdo_calibration(do_pdo_calibration, pdoCalibFile)

    # get the run number from input file
    run_number = utils.run_number_from_file(input_file)
    vmmana.run_number = run_number

    vmmana.load_tree("vmm2", input_file)

    # determine how many events to process
    n_total = vmmana.tree.GetEntries()
    n_to_process = n_total
    if nevents >= 0 and not (nevents>n_total) :
        n_to_process = nevents

    for ievent, event in enumerate(vmmana.tree) :
        if ievent%5000==0 :
            info("Processing entry %d/%d"%(ievent,n_to_process))
        vmmana.clear_containers()
        process_event(vmmana, ievent, event)

