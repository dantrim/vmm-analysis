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

    
