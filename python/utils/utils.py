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

def get_tchamber_strip(vmm_id, vmm_channel) :
    if (vmm_id<0) or (vmm_channel<0) :
        print "get_tchamber_strip    invalid vmm id and/or channel"
        sys.exit()
    
    strip = 0
    if vmm_id%2==0 :
        if vmm_channel== 0  :  strip= 32
        if vmm_channel== 1  :  strip= 33
        if vmm_channel== 2  :  strip= 31
        if vmm_channel== 3  :  strip= 34
        if vmm_channel== 4  :  strip= 30
        if vmm_channel== 5  :  strip= 35
        if vmm_channel== 6  :  strip= 29
        if vmm_channel== 7  :  strip= 36
        if vmm_channel== 8  :  strip= 28
        if vmm_channel== 9  :  strip= 37
        if vmm_channel== 10 :  strip= 27
        if vmm_channel== 11 :  strip= 38
        if vmm_channel== 12 :  strip= 26
        if vmm_channel== 13 :  strip= 39
        if vmm_channel== 14 :  strip= 25
        if vmm_channel== 15 :  strip= 40
        if vmm_channel== 16 :  strip= 24
        if vmm_channel== 17 :  strip= 41
        if vmm_channel== 18 :  strip= 23
        if vmm_channel== 19 :  strip= 42
        if vmm_channel== 20 :  strip= 22
        if vmm_channel== 21 :  strip= 43
        if vmm_channel== 22 :  strip= 21
        if vmm_channel== 23 :  strip= 44
        if vmm_channel== 24 :  strip= 20
        if vmm_channel== 25 :  strip= 45
        if vmm_channel== 26 :  strip= 19
        if vmm_channel== 27 :  strip= 46
        if vmm_channel== 28 :  strip= 18
        if vmm_channel== 29 :  strip= 47
        if vmm_channel== 30 :  strip= 17
        if vmm_channel== 31 :  strip= 48
        if vmm_channel== 32 :  strip= 16
        if vmm_channel== 33 :  strip= 49
        if vmm_channel== 34 :  strip= 15
        if vmm_channel== 35 :  strip= 50
        if vmm_channel== 36 :  strip= 14
        if vmm_channel== 37 :  strip= 51
        if vmm_channel== 38 :  strip= 13
        if vmm_channel== 39 :  strip= 52
        if vmm_channel== 40 :  strip= 12
        if vmm_channel== 41 :  strip= 53
        if vmm_channel== 42 :  strip= 11
        if vmm_channel== 43 :  strip= 54
        if vmm_channel== 44 :  strip= 10
        if vmm_channel== 45 :  strip= 55
        if vmm_channel== 46 :  strip= 9 
        if vmm_channel== 47 :  strip= 56
        if vmm_channel== 48 :  strip= 8 
        if vmm_channel== 49 :  strip= 57
        if vmm_channel== 50 :  strip= 7 
        if vmm_channel== 51 :  strip= 58
        if vmm_channel== 52 :  strip= 6 
        if vmm_channel== 53 :  strip= 59
        if vmm_channel== 54 :  strip= 5 
        if vmm_channel== 55 :  strip= 60
        if vmm_channel== 56 :  strip= 4 
        if vmm_channel== 57 :  strip= 61
        if vmm_channel== 58 :  strip= 3 
        if vmm_channel== 59 :  strip= 62
        if vmm_channel== 60 :  strip= 2 
        if vmm_channel== 61 :  strip= 63
        if vmm_channel== 62 :  strip= 1 
        if vmm_channel== 63 :  strip= 64
    else :
        if vmm_channel== 0  :   strip=   128  
        if vmm_channel== 1  :   strip=   65   
        if vmm_channel== 2  :   strip=   127  
        if vmm_channel== 3  :   strip=   66   
        if vmm_channel== 4  :   strip=   126  
        if vmm_channel== 5  :   strip=   67   
        if vmm_channel== 6  :   strip=   125  
        if vmm_channel== 7  :   strip=   68   
        if vmm_channel== 8  :   strip=   124  
        if vmm_channel== 9  :   strip=   69   
        if vmm_channel== 10 :   strip=   123  
        if vmm_channel== 11 :   strip=   70   
        if vmm_channel== 12 :   strip=   122  
        if vmm_channel== 13 :   strip=   71   
        if vmm_channel== 14 :   strip=   121  
        if vmm_channel== 15 :   strip=   72   
        if vmm_channel== 16 :   strip=   120  
        if vmm_channel== 17 :   strip=   73   
        if vmm_channel== 18 :   strip=   119  
        if vmm_channel== 19 :   strip=   74   
        if vmm_channel== 20 :   strip=   118  
        if vmm_channel== 21 :   strip=   75   
        if vmm_channel== 22 :   strip=   117  
        if vmm_channel== 23 :   strip=   76   
        if vmm_channel== 24 :   strip=   116  
        if vmm_channel== 25 :   strip=   77   
        if vmm_channel== 26 :   strip=   115  
        if vmm_channel== 27 :   strip=   78   
        if vmm_channel== 28 :   strip=   114  
        if vmm_channel== 29 :   strip=   79   
        if vmm_channel== 30 :   strip=   113  
        if vmm_channel== 31 :   strip=   80   
        if vmm_channel== 32 :   strip=   112  
        if vmm_channel== 33 :   strip=   81   
        if vmm_channel== 34 :   strip=   111  
        if vmm_channel== 35 :   strip=   82   
        if vmm_channel== 36 :   strip=   110  
        if vmm_channel== 37 :   strip=   83   
        if vmm_channel== 38 :   strip=   109  
        if vmm_channel== 39 :   strip=   84   
        if vmm_channel== 40 :   strip=   108  
        if vmm_channel== 41 :   strip=   85   
        if vmm_channel== 42 :   strip=   107  
        if vmm_channel== 43 :   strip=   86   
        if vmm_channel== 44 :   strip=   106  
        if vmm_channel== 45 :   strip=   87   
        if vmm_channel== 46 :   strip=   105  
        if vmm_channel== 47 :   strip=   88   
        if vmm_channel== 48 :   strip=   104  
        if vmm_channel== 49 :   strip=   89   
        if vmm_channel== 50 :   strip=   103  
        if vmm_channel== 51 :   strip=   90   
        if vmm_channel== 52 :   strip=   102  
        if vmm_channel== 53 :   strip=   91   
        if vmm_channel== 54 :   strip=   101  
        if vmm_channel== 55 :   strip=   92   
        if vmm_channel== 56 :   strip=   100  
        if vmm_channel== 57 :   strip=   93   
        if vmm_channel== 58 :   strip=   9    
        if vmm_channel== 59 :   strip=   94   
        if vmm_channel== 60 :   strip=   8    
        if vmm_channel== 61 :   strip=   95   
        if vmm_channel== 62 :   strip=   7    
        if vmm_channel== 63 :   strip=   96   

    if vmm_id%2==0 and (vmm_id==2 or vmm_id==6) :
        strip = strip + 128
    elif (not (vmm_id%2==0)) and (vmm_id==3 or vmm_id==7) :
        strip = strip + 128

    return strip

   # if(vmm_id%2==0 && (vmm_id==2 || vmm_id==3)) {
   #     strip+=128;
   # }
   # else if(!(vmm_id%2==0) && (vmm_id==6 || vmm_id==7)) {
   #     strip+=128;
   # }
