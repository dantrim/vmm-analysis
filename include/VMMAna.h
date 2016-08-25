//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 23 13:25:28 2016 by ROOT version 6.04/00
// from TTree vmm2/vmm2
// found on file: run_0008.root
//////////////////////////////////////////////////////////

#ifndef VMMAna_h
#define VMMAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// me
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "hit.h"

// Header file for the classes stored in the TTree if any.
#include "vector"

class VMMAna : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain


    // MY CLASSES
    vector<Hit*> hits;
    void build_hits();

    // functions
    void init_canvases();
    void init_raw_histos();
    void init_hit_map_histos();
    void init_duplicate_histos();



    void loadContainers();
    void clearContainers();
    void fill_raw_quantities();
    void performDuplicateRemoval();

    void sort_by_bcid();

    //// counters
    bool duplicate_found;
    int n_duplicate_per_channel;
    int n_duplicate_per_event;
    int n_crossfire_per_event;
    int n_cluster_before_duplicate_removal;


    void draw_raw_histos();

    void draw_duplicate_histos();

    int TChamberMapping(int vmm_id, int vmm_chan);

    // event dudes
    vector<int> _chip;
    vector< vector<int> > _pdo;
    vector< vector<int> > _tdo;
    vector< vector<int> > _bcid;
    vector< vector<int> > _gray;
    vector< vector<int> > _channel;
    vector< vector<int> > _threshold;

    // clusters
    //tmp
    vector<int> tmp_cluster_channel;
    vector<int> tmp_cluster_strips;
    vector<double> tmp_cluster_pdos;
    vector<int> tmp_cluster_bcids;
    vector<int> tmp_cluster_grays;
    vector<int> tmp_cluster_chips;

    vector<int> cluster_channel;
    vector<int> cluster_strips;
    vector<double> cluster_pdos;
    vector<int> cluster_bcids;
    vector<int> cluster_grays;
    vector<int> cluster_chips;

    vector<vector<int> > event_cluster_channel;
    vector<vector<int> > event_cluster_strips;
    vector<vector<double> > event_cluster_pdos;
    vector<vector<int> > event_cluster_bcids;
    vector<vector<int> > event_cluster_grays;
    vector<vector<int> > event_cluster_chips;

    
    


    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    // PLOTTING THINGS [BEGIN]
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    // canvas
    TCanvas* c_charge; // canvas for channel/cluster charges
    TCanvas* c_hit_maps; // canvas for showing the VMM hit maps
    TCanvas* c_chamber_hits; // canvas for showing the chamber raw hit maps
    TCanvas* c_event_hits; // canvas showing the 2D strip vs pdo hit map
    TCanvas* c_duplicate_channel; // canvas for showing duplicate channels per VMM
    TCanvas* c_duplicate_strip;
    TCanvas* c_total_duplicates; // canvas for showing total # duplicates per VMM per event


    ////////////////// 1D HISTOS ////////////////////

    // charge histograms
    vector<TH1F*> h_pdo_raw; // charge distribution, all channels, no calib

    // hit map histograms
    vector<TH1F*> h_vmm_hits;

    // chamber raw hit maps
    vector<TH1F*> h_chamber_hits;

    // duplicates per channel
    vector<TH1F*> h_duplicate_chan;
    vector<TH1F*> h_duplicate_strip;
    TH1F* h_total_duplicates;
    TH1F* h_total_cross_duplicates;

    ////////////////// 2D HISTOS ////////////////////
    vector<TH2F*> h2_event_hit_map; // strip vs pdo

    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    // PLOTTING THINGS [END]
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventFAFA;
   vector<int>     *triggerTimeStamp;
   vector<int>     *triggerCounter;
   vector<int>     *chip;
   vector<int>     *eventSize;
   vector<vector<int> > *tdo;
   vector<vector<int> > *pdo;
   vector<vector<int> > *flag;
   vector<vector<int> > *threshold;
   vector<vector<int> > *bcid;
   vector<vector<int> > *grayDecoded;
   vector<vector<int> > *channel;

   // List of branches
   TBranch        *b_eventFAFA;   //!
   TBranch        *b_triggerTimeStamp;   //!
   TBranch        *b_triggerCounter;   //!
   TBranch        *b_chip;   //!
   TBranch        *b_eventSize;   //!
   TBranch        *b_tdo;   //!
   TBranch        *b_pdo;   //!
   TBranch        *b_flag;   //!
   TBranch        *b_threshold;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_grayDecoded;   //!
   TBranch        *b_channel;   //!

   VMMAna(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~VMMAna() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(VMMAna,0);
};

#endif

#ifdef VMMAna_cxx
//void VMMAna::Init(TTree *tree)
//{
//   // The Init() function is called when the selector needs to initialize
//   // a new tree or chain. Typically here the branch addresses and branch
//   // pointers of the tree will be set.
//   // It is normally not necessary to make changes to the generated
//   // code, but the routine can be extended by the user if needed.
//   // Init() will be called many times when running on PROOF
//   // (once per file to be processed).
//
//   // Set object pointer
//   triggerTimeStamp = 0;
//   triggerCounter = 0;
//   chip = 0;
//   eventSize = 0;
//   tdo = 0;
//   pdo = 0;
//   flag = 0;
//   threshold = 0;
//   bcid = 0;
//   grayDecoded = 0;
//   channel = 0;
//   // Set branch addresses and branch pointers
//   if (!tree) return;
//   fChain = tree;
//   fChain->SetMakeClass(1);
//
//   fChain->SetBranchAddress("eventFAFA", &eventFAFA, &b_eventFAFA);
//   fChain->SetBranchAddress("triggerTimeStamp", &triggerTimeStamp, &b_triggerTimeStamp);
//   fChain->SetBranchAddress("triggerCounter", &triggerCounter, &b_triggerCounter);
//   fChain->SetBranchAddress("chip", &chip, &b_chip);
//   fChain->SetBranchAddress("eventSize", &eventSize, &b_eventSize);
//   fChain->SetBranchAddress("tdo", &tdo, &b_tdo);
//   fChain->SetBranchAddress("pdo", &pdo, &b_pdo);
//   fChain->SetBranchAddress("flag", &flag, &b_flag);
//   fChain->SetBranchAddress("threshold", &threshold, &b_threshold);
//   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
//   fChain->SetBranchAddress("grayDecoded", &grayDecoded, &b_grayDecoded);
//   fChain->SetBranchAddress("channel", &channel, &b_channel);
//}

Bool_t VMMAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef VMMAna_cxx
