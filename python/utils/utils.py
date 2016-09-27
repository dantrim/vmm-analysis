#!/usr/bin/env python

import os

def file_ok(input_file="") :
    if input_file=="" :
        print "utils::file_ok    Provided input file is \"\""
        return False

    if os.path.exists(input_file) :
        print "utils::file_ok    Found provided file: %s"%input_file
        return True
    else :
        print "utils::file_ok    Unable to find provided file: %s"%input_file
        return False

def run_number_from_file(input_file) :
    """
    returns the run number as a string.
    expects file format of "run_XXXX.root"
    """
    run = input_file.strip()
    run = run.split(".root")[0]
    run = run.split("_")[1]
    return run

    
