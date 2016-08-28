#include "tools.hh"
#include <cstdio>
#include <cstdlib>

//std/stl
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

//ROOT
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"

//vmm
#include "t6t7_ana.hh"

void help()
{
    cout << "Options:" << endl
         << "   -i|--input          name of input file" << endl
         << "   -n|--nEvents        number of events to process (default: -1, all events)" << endl
         << "   -d|--dbg            run in debug/verbose mode" << endl
         << "   -h|--help           print this help message" << endl
    << endl;
}


int main(int argc, char* argv[])
{
    gSystem->Load("/Users/dantrim/workarea/NSW/myanalysis/vmm-analysis/vector_lib/libMylib");

    int nEvt            = -1;
    int dbg             = 0;
    string inputFile    = "";
    string calibFile    = "";

    cout << " -------------------------------- " << endl;
    cout << "            T6/T7 ana             " << endl;
    cout << " -------------------------------- " << endl;

    int optin(1);
    while (optin < argc) {
        string opt = argv[optin];
        if      (opt=="-i" || opt=="--input")   { inputFile = argv[++optin]; }
        else if (opt=="-d" || opt=="--dbg")     { dbg = true; }
        else if (opt=="-n" || opt=="--nEvents") { nEvt = atoi(argv[++optin]); }
        else if (opt=="-c" || opt=="--calibFile") { calibFile = argv[++optin]; }
        else if (opt=="-h" || opt=="--help")    { help(); return 0; }
        else {
            cout << "Unknown command line argument : '" << opt << "'" << endl;
            help();
            return 1;
        }
        optin++;
    }

    bool do_calib = false;

    if(inputFile=="") {
        cout << "You have not provided an input file. Exitting." << endl;
        help();
        return 1;
    }
    bool exists = std::ifstream(inputFile).good();
    if(!exists) {
        cout << "Problem locating input file. Exitting." << endl;
        return 1;
    }

    if(calibFile!="") do_calib = true;

    int run_number = -1;
    size_t find = inputFile.find("run_");
    string run = inputFile.substr(find+4, find+4);
    run_number = stoi(run);

    TChain* chain = new TChain("vmm2");
    TFile* file = new TFile(inputFile.c_str());
    chain = static_cast<TChain*>(file->Get("vmm2"));
    Long64_t nEntries = chain->GetEntries(); 
    chain->ls();

    t6t7::VMMAna* vmm_ana = new t6t7::VMMAna();
    vmm_ana->setDebug(dbg);
    vmm_ana->setChain(chain);
    vmm_ana->setRunNumber(run_number);
    vmm_ana->setDoCalib(do_calib);
    vmm_ana->setCalibFile(calibFile);
    

    if(nEvt<0) nEvt = nEntries;
    cout << endl;
    cout << " Total entries      : " << nEntries << endl;
    cout << " Entries to process : " << nEvt << endl;

    chain->Process(vmm_ana, "", nEvt, 0);


    cout << " ------------ VMMAna job done --------------- " << endl;
    delete chain;
    delete vmm_ana;


//  int in = 0;
//  if (argc > 1) in = atoi(argv[1]);
//  int testint = testfunc(in);
//  printf("bonjour et %i\n", testint);
  return 0;
}
