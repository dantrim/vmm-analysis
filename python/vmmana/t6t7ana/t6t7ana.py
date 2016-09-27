#!/usr/bin/env python

from optparse import OptionParser
import os

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)
#r.gStyle.SetOptStat(False)

import sys
sys.path.append("../../utils/")

# vmmana
import cluster as vmmcluster
import hit as vmmhit
import utils
import plot_utils as pu

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
        self.tdo_constants = {}
        self.mean_tdo_constant = 0 # for when we don't have constants for ALL channels

        self.do_pdo_calib = False
        self.pdo_file = ""
        self.pdo_constants = {}

        self.tree_file = None
        self.tree = None

        # container for all hits
        self.hits = []
    
        # container for all clusters
        self.clusters = []

        # canvases with plots
        self.canvases = []

        # histograms
        self.histograms = {}

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

    def get_tdo_constant(self, chip, channel) :
        for loaded_chip in self.tdo_constants.keys() :
            if loaded_chip == chip :
                return self.tdo_constants[loaded_chip][channel]
        return self.mean_tdo_constant
        

    def load_tdo_calibration(self, do_calib=False, calib_file="") :
        """
        loads in the tdo calibration constants from *.txt file
        with 3 columns per line:
            [1] chip number
            [2] vmm-channel
            [3] calibration constant [ns/ADC cnts]
        """
        self.do_tdo_calib = do_calib
        if do_calib :
            self.tdo_file = calib_file
            lines = open(calib_file).readlines()
            
            mean_constant = 0.0
            n_constants = 0
            chips_in_file = []
            for line in lines :
                if not line : continue
                line = line.strip()
                if line.startswith("#") : continue
                line = line.split()
                chip = line[0]
                channel = line[1]
                constant = line[2]
                mean_constant += float(constant)
                n_constants += 1
                if chip not in chips_in_file :
                    chips_in_file.append(chip)

            for chip in chips_in_file :
                self.tdo_constants[int(chip)] = {}

            for line in lines :
                if not line : continue
                if line.startswith("#") : continue
                line = line.split()
                chip = line[0]
                channel = line[1]
                constant = line[2]
                self.tdo_constants[int(chip)][int(channel)] = float(constant)

            mean_constant = mean_constant / n_constants
            self.mean_tdo_constant = mean_constant
    

    def load_pdo_calibration(self, do_calib=False, calib_file="") :
        print "load_pdo_calibration    NOT YET IMPLEMENTED"
        sys.exit()

    def add_canvas(self, c) :
        if c==None :
            info("add_canvas    Input canvas is Null!")
            sys.exit()
        self.canvases.append(c)

def load_chamber_hits(vmmdata) :

    out_hits = []
    #for chip in vmmdata.chip :
    n_chip = len(vmmdata.chip)
    for ichip in xrange(n_chip) :
        n_hits = len(vmmdata.tdo[ichip])
        for ihit in xrange(n_hits) :
            h = vmmhit.Hit()
            chip_number = vmmdata.chip[ichip]
            channel = vmmdata.channel[ichip][ihit]
            bcid = vmmdata.bcid[ichip][ihit]
            gray = vmmdata.grayDecoded[ichip][ihit]
            strip = utils.get_tchamber_strip(chip_number, channel)

            # for T7, rotate the strips
            if chip_number>=4 and chip_number<8 :
                strip = 258 - strip 

            threshold = vmmdata.threshold[ichip][ihit]
            pdo = vmmdata.pdo[ichip][ihit]
            tdo = vmmdata.tdo[ichip][ihit]

            chamber = 0
            if chip_number>=0 and chip_number<4 :
                chamber = 0 
            elif chip_number>=4 and chip_number<8 :
                chamber = 1

            h.setChamber(chamber)
            h.fill(chip_number, bcid, gray, channel, strip, threshold, pdo, tdo)
            out_hits.append(h)
    return out_hits

def remove_duplicate_hits(vmm) :


    hits_to_remove = []
    for i, hiti in enumerate(vmm.hits) :
        for j, hitj in enumerate(vmm.hits) :
            if j>i :
                if hiti==hitj :
                    hits_to_remove.append(hiti)
                    break # removed i, no longer compare it to j

    new_hits = []
    for hit in vmm.hits :
        add_hit = True
        for duplicate_hit in hits_to_remove :
            if hit==duplicate_hit :
                add_hit = False
                break
        if add_hit :
            new_hits.append(hit)

    vmm.hits = []
    vmm.hits = new_hits

###############################################################################
# event loop
def process_event(vmm, event_no, vmmdata) :

    vmm.hits = load_chamber_hits(vmmdata)
    # sort hits earliest to latest in bcid
    vmm.hits.sort(key=lambda x : x.gray, reverse=False) 
    
    fill_raw_histograms(vmm)

    remove_duplicate_hits(vmm)
    vmm.hits.sort(key=lambda x : x.gray, reverse=False)



###############################################################################
# histograms
def initialize_histograms(vmm) :

    # clear the histograms
    vmm.histograms = {}

    ################################################
    # initialize raw charge histograms
    raw_pdo_per_chip = []
    x_title = "PDO [ADC counts]"
    y_title = "# of entries"
    for ichip in xrange(8) :
        h = r.TH1F("h_raw_pdo_chip%d"%ichip, "VMM %d Charge Distribution (all chan, no calib)"%ichip, 400, 0, 1200)
        h.GetXaxis().SetTitle(x_title)
        h.GetYaxis().SetTitle(y_title)
        raw_pdo_per_chip.append(h)
    vmm.histograms["raw_pdo_per_chip"] = raw_pdo_per_chip

    ################################################
    # initialize tdo histograms
    raw_tdo_per_chip = []
    x_title = "TDO [ADC counts]"
    y_title = "# of entries"
    for ichip in xrange(8) :
        h = r.TH1F("h_raw_tdo_chip%d"%ichip, "VMM %d TDO Distribution (all, chan, no calib)"%ichip, 250, 0, 800)
        h.GetXaxis().SetTitle(x_title)
        h.GetYaxis().SetTitle(y_title)
        raw_tdo_per_chip.append(h)
    vmm.histograms["raw_tdo_per_chip"] = raw_tdo_per_chip
    
    ################################################
    # initialize channel hit maps
    vmm_channel_hits = []
    for ichip in xrange(8) :
        h = r.TH1F("h_vmm_channel_hits_chip%d"%ichip,"",63,0,63)
        h.SetTitle("VMM %d Raw Channel Hits"%ichip)
        x_title = "VMM %d channel"%ichip
        y_title = "# of entries"
        h.GetXaxis().SetTitle(x_title)
        h.GetYaxis().SetTitle(y_title)
        vmm_channel_hits.append(h)
    vmm.histograms["vmm_channel_hits"] = vmm_channel_hits

    ################################################
    # initialize strip hits
    raw_strip_hits = []
    chamber_names = ["T6", "T7"]
    x_title = "MicroMegas strip"
    y_title = "# of entries"
    for ichamber, chamber in enumerate(chamber_names) :
        h = r.TH1F("h_raw_strip_hits_chamber%s"%chamber, "", 256, 1, 257)
        h.GetXaxis().SetTitle(x_title)
        h.GetYaxis().SetTitle(y_title)
        raw_strip_hits.append(h)
    vmm.histograms["raw_strip_hits"] = raw_strip_hits

    #################################################
    # raw time histograms
    if vmm.do_tdo_calib :
        raw_hit_times_per_chip = []
        for ichip in xrange(8) :
            h = r.TH1F("h_raw_times_chip%d"%ichip, "", 200,0,200)
            h.SetTitle("VMM %d Raw Hit Times"%ichip)
            x_title = "Time [ns]"
            y_title = "# of entries"
            h.GetXaxis().SetTitle(x_title)
            h.GetYaxis().SetTitle(y_title)
            raw_hit_times_per_chip.append(h)
        vmm.histograms["raw_hit_times_chip"] = raw_hit_times_per_chip


        raw_hit_times_per_chamber = []
        chamber_names = ["T6", "T7"]
        x_title = "Time [ns]"
        y_title = "# of entries"
        for ichanber, chamber in enumerate(chamber_names) :
            h = r.TH1F("h_raw_times_chamber%s"%chamber, "", 200, 0, 200)
            h.GetXaxis().SetTitle(x_title)
            h.GetYaxis().SetTitle(y_title)
            raw_hit_times_per_chamber.append(h)
        vmm.histograms["raw_hit_times_chamber"] = raw_hit_times_per_chamber
         

def fill_raw_histograms(vmm) :

    for hit in vmm.hits :
        if hit.chip_number>=8 or hit.chip_number<0 : continue
        if not (hit.bcid>0) : continue
        if hit.pdo>0 and hit.pdo<1024 and hit.threshold==1 :
            vmm.histograms["raw_pdo_per_chip"][hit.chip_number].Fill(hit.pdo) 
            vmm.histograms["raw_tdo_per_chip"][hit.chip_number].Fill(hit.tdo)
            vmm.histograms["vmm_channel_hits"][hit.chip_number].Fill(hit.channel)
            vmm.histograms["raw_strip_hits"][hit.chamber_number].Fill(hit.strip)
            if vmm.do_tdo_calib :
                tac = 125 #ns
                raw_hit_time = hit.tdo * vmm.get_tdo_constant(hit.chip_number, hit.channel)
                vmm.histograms["raw_hit_times_chip"][hit.chip_number].Fill(raw_hit_time)
                vmm.histograms["raw_hit_times_chamber"][hit.chamber_number].Fill(raw_hit_time)


def draw_histograms(vmm) :
    draw_raw_histograms(vmm)


def draw_raw_histograms(vmm) :

    ###############################################
    # raw charge histograms
    key = "raw_pdo_per_chip"
    histos = vmm.histograms[key]
    name = "%s_%s"%(key, vmm.run_number)
    c = pu.VMMCanvas("%s"%name, 1400, 800) 
    c.canvas().Divide(4,2)
    for i, histo in enumerate(histos) :
        c.canvas().cd(i+1)
        histo.Draw()
    vmm.add_canvas(c)

    ###############################################
    # raw tdo histograms
    key = "raw_tdo_per_chip"
    histos = vmm.histograms[key]
    name = "%s_%s"%(key, vmm.run_number)
    c = pu.VMMCanvas("%s"%name, 1400, 800)
    c.canvas().Divide(4,2)
    for i, histo in enumerate(histos) :
        c.canvas().cd(i+1)
        histo.Draw()
    vmm.add_canvas(c) 

    ###############################################
    # vmm channel hits 
    key = "vmm_channel_hits"
    histos = vmm.histograms[key]
    name = "%s_%s"%(key, vmm.run_number)
    c = pu.VMMCanvas("%s"%name, 1400, 800) 
    c.canvas().Divide(4,2)
    for i, histo in enumerate(histos) :
        c.canvas().cd(i+1)
        histo.Draw()
    vmm.add_canvas(c)

    ###############################################
    # chamber strip hits
    key = "raw_strip_hits"
    histos = vmm.histograms[key]
    name = "%s_%s"%(key, vmm.run_number)
    c = pu.VMMCanvas("%s"%name, 1400, 800) 
    c.canvas().Divide(1,2)
    for i, histo in enumerate(histos) :
        c.canvas().cd(i+1)
        histo.Draw()
    vmm.add_canvas(c)

    ###############################################
    if vmm.do_tdo_calib :

        # chip hit times raw
        key = "raw_hit_times_chip"
        histos = vmm.histograms[key]
        name = "%s_%s"%(key, vmm.run_number)
        c = pu.VMMCanvas("%s"%name, 1400, 800)
        c.canvas().Divide(4,2)
        for i, histo in enumerate(histos) :
            c.canvas().cd(i+1)
            histo.Draw()
        vmm.add_canvas(c)

        # chamber hit times raw
        key = "raw_hit_times_chamber"
        histos = vmm.histograms[key]
        name = "%s_%s"%(key, vmm.run_number)
        c = pu.VMMCanvas("%s"%name, 800, 600)
        c.canvas().Divide(1,2)
        for i, histo in enumerate(histos) :
            c.canvas().cd(i+1)
            histo.Draw()
        vmm.add_canvas(c) 

def write_plots(vmm, output_dir = "./") :
    for can in vmm.canvases :
        if not output_dir.endswith("/") : output_dir = output_dir + "/"
        save_name = "%s%s.eps"%(output_dir, can.name)
        can.canvas().SaveAs(save_name)

###############################################################################
if __name__=="__main__" :
    print "starting..."

    parser = OptionParser()
    parser.add_option("-i", "--input", default="")
    parser.add_option("-d", "--debug", action="store_true", default=False)
    parser.add_option("-n", "--nEvents", default=-1)
    parser.add_option("--tdoCalibFile", default="")
    parser.add_option("--pdoCalibFile", default="")
    parser.add_option("-o", "--outDir", default="./")
    (options, args) = parser.parse_args()

    input_file = options.input 
    dbg = options.debug
    nevents = int(options.nEvents)
    tdoCalibFile = options.tdoCalibFile
    pdoCalibFile = options.pdoCalibFile
    output_dir = options.outDir

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
    vmmana.load_tdo_calibration(do_tdo_calibration, tdoCalibFile)
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

    # initialize histograms
    initialize_histograms(vmmana)

    for ievent, event in enumerate(vmmana.tree) :
        if (ievent % 5000) == 0 :
            info("Processing entry %d/%d"%(ievent,n_to_process))
        if ievent >= n_to_process :
            break
        vmmana.clear_containers()
        process_event(vmmana, ievent, event)

    # draw_histograms
    draw_histograms(vmmana)

    # save the plots
    info("Storing plots in directory: %s"%output_dir)
    write_plots(vmmana, output_dir)

