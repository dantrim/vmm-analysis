#define VMMAna_cxx
// The class definition in VMMAna.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("VMMAna.C")
// root> T->Process("VMMAna.C","some options")
// root> T->Process("VMMAna.C+")
//

#include "VMMAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;


void VMMAna::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}
void VMMAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   triggerTimeStamp = 0;
   triggerCounter = 0;
   chip = 0;
   eventSize = 0;
   tdo = 0;
   pdo = 0;
   flag = 0;
   threshold = 0;
   bcid = 0;
   grayDecoded = 0;
   channel = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventFAFA", &eventFAFA, &b_eventFAFA);
   fChain->SetBranchAddress("triggerTimeStamp", &triggerTimeStamp, &b_triggerTimeStamp);
   fChain->SetBranchAddress("triggerCounter", &triggerCounter, &b_triggerCounter);
   fChain->SetBranchAddress("chip", &chip, &b_chip);
   fChain->SetBranchAddress("eventSize", &eventSize, &b_eventSize);
   fChain->SetBranchAddress("tdo", &tdo, &b_tdo);
   fChain->SetBranchAddress("pdo", &pdo, &b_pdo);
   fChain->SetBranchAddress("flag", &flag, &b_flag);
   fChain->SetBranchAddress("threshold", &threshold, &b_threshold);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("grayDecoded", &grayDecoded, &b_grayDecoded);
   fChain->SetBranchAddress("channel", &channel, &b_channel);

    /////////////////////////////////////////////////////////////////////

    // initialize the raw histograms
    init_raw_histos();

    // initialize hit map histograms
    init_hit_map_histos();

    // initialize duplicate histos
    init_duplicate_histos();
}

void VMMAna::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t VMMAna::Process(Long64_t entry)
{
    loadContainers();

    fChain->GetEntry(entry);

    cout << " *** Processing entry " << entry << " *** " << endl;

    build_hits();

    fill_raw_quantities();

    sort_by_bcid();

    duplicate_found = false;
    n_duplicate_per_channel = 0;
    n_duplicate_per_event = 0;
    n_crossfire_per_event = 0;
    n_cluster_before_duplicate_removal = tmp_cluster_strips.size();

    performDuplicateRemoval();


    clearContainers();
    return kTRUE;
}

void VMMAna::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
    // draw the histograms

}

void VMMAna::Terminate()
{

    cout << "VMMAna::Terminate" << endl;
    draw_raw_histos();
    draw_duplicate_histos();


}

//------------ USER DEFINED ----------------- //
void VMMAna::loadContainers()
{
    _chip = *chip;
    _pdo = *pdo;
    _tdo = *tdo;
    _bcid = *bcid;
    _gray = *grayDecoded;
    _channel = *channel;
    _threshold = *threshold;

    return;
}
void VMMAna::clearContainers()
{
    _chip.clear();
    _pdo.clear();
    _bcid.clear();
    _gray.clear();
    _channel.clear();
    _threshold.clear();

    tmp_cluster_strips.clear();
    tmp_cluster_pdos.clear();
    tmp_cluster_bcids.clear();
    tmp_cluster_grays.clear();
    tmp_cluster_chips.clear();

    cluster_strips.clear();
    cluster_pdos.clear();
    cluster_bcids.clear();
    cluster_grays.clear();
    cluster_chips.clear();

    event_cluster_strips.clear();
    event_cluster_pdos.clear();
    event_cluster_bcids.clear();
    event_cluster_grays.clear();
    event_cluster_chips.clear();

    
}
void VMMAna::init_canvases()
{
    // charge histograms

}


void VMMAna::init_raw_histos()
{
    // histogram of raw charge over all channels
    stringstream name;
    for(int i = 0; i < 8; i++) {
        name << "VMM " << i << " Charge Distribution (all chan, no calib)";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 400, 0, 1200);
        h->GetXaxis()->SetTitle("PDO [ADC counts]");
        h->GetYaxis()->SetTitle("# of Entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_pdo_raw.push_back(h);
        name.str("");
    }

    //c_event_hits = new TCanvas("c_event_hits", "", 1600, 800);
    //c_event_hits->Divide(1,2);
}

void VMMAna::init_hit_map_histos()
{
    stringstream name;
    for(int i = 0; i < 8; i++) {
        name << "h_vmm_" << i << "_hits";
        TH1F* h = new TH1F(name.str().c_str(), "", 63, 0, 63);
        name.str("");
        name << "chip " << i << " - VMM channel";
        h->GetXaxis()->SetTitle(name.str().c_str());
        name.str("");
        name << "# of Entries";
        h->GetYaxis()->SetTitle(name.str().c_str());
        h->GetYaxis()->SetTitleOffset(1.5);
        
        h_vmm_hits.push_back(h);
    }

    name.str("");
    for(int i = 0; i < 2; i++) {
        name << "T";
        if(i==0) name << "6";
        else { name << "7"; }
        name << " Chamber Raw Hit Map";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 256,0,256);
        h->GetXaxis()->SetTitle("MicroMegas Strip");
        h->GetYaxis()->SetTitle("# of Entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_chamber_hits.push_back(h);
        name.str("");
    }

    name.str("");
    for(int i = 0; i < 2; i++) {
        name << "T";
        if(i==0) name << "6";
        else { name << "7"; }
        name << " Event Hit Map";
        TH2F* h = new TH2F(name.str().c_str(), name.str().c_str(), 256, 0, 256, 1100, 0, 1100);
        h->SetMarkerStyle(26);
        h->SetMarkerSize(1.2);
        h->GetXaxis()->SetTitle("Strip");
        h->GetYaxis()->SetTitle("PDO [ADC counts]");
        h2_event_hit_map.push_back(h);
        name.str("");
    }
}

void VMMAna::init_duplicate_histos()
{
    stringstream name;
    for(int i = 0; i < 8; i++) {
        name << "# of Duplicate Entries/Channel (VMM-" << i << ")";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 64,0,64);
        h->GetXaxis()->SetTitle("VMM channel");
        h->GetYaxis()->SetTitle("# of Duplicates");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_duplicate_chan.push_back(h);
        name.str("");
    }

    name.str("");
    for(int i = 0; i < 2; i++) {
        name << "# of Duplicates per T" << ( i==0 ? "6" : "7" ) << " Chamber";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 256, 0, 256);
        h->GetXaxis()->SetTitle("T-Chamber Strip");
        h->GetYaxis()->SetTitle("# of Duplicates");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_duplicate_strip.push_back(h);
        name.str("");
    }

    name.str("");
    name << "Total # of Duplicate Entires per VMM Channel";
    h_total_duplicates = new TH1F(name.str().c_str(), name.str().c_str(), 20, 0, 20);
    h_total_duplicates->GetXaxis()->SetTitle("# of Duplicates");
    h_total_duplicates->GetYaxis()->SetTitle("# of Entries");
    h_total_duplicates->GetYaxis()->SetTitleOffset(1.5);

    name.str("");
    name << "Total # of Cross-Channel Duplicates per Event";
    h_total_cross_duplicates = new TH1F(name.str().c_str(), name.str().c_str(), 20, 0, 20);
    h_total_duplicates->GetXaxis()->SetTitle("# of Duplicates");
    h_total_duplicates->GetYaxis()->SetTitle("# of Entries");
    h_total_duplicates->GetYaxis()->SetTitleOffset(1.5);

}

void VMMAna::build_hits()
{
    hits.clear();
    Hit* hit = new Hit();

  //  for(int ichip = 0; ichip < (int)_chip.size(); ichip++) {
  //      for(int ichan = 0; ichan < (int)_pdo.at(ichip).size(); ichan++) {
  //          Hit* hit = new Hit();
  //          int chip_ = _chip[ichip];
  //          int pdo_ = _pdo[ichip][ichan];
  //          int tdo_ = _tdo[ichip][ichan];
  //          int bcid_ = _bcid[ichip][ichan];
  //          int gray_ = _gray[ichip][ichan];
  //          int channel_ = _channel[ichip][ichan];
  //          int threshold_ = _threshold[ichip][ichan];

  //          int strip_ = TChamberMapping(chip_, channel_);

  //          hit->fill(chip_,
  //                      bcid_,
  //                      gray_,
  //                      channel_,
  //                      strip_,
  //                      threshold_,
  //                      pdo_,
  //                      tdo_);
  //          hit->print();
  //          hits.push_back(hit);
  //      }
  //  }

}

void VMMAna::fill_raw_quantities()
{

    h2_event_hit_map.at(0)->Reset();
    h2_event_hit_map.at(1)->Reset();
    bool has_hits = false;

    for(int i = 0; i < (int)_pdo.size(); i++) {
        if(_chip[i]>=8) continue;
        for(int j = 0; j < (int)_pdo[i].size(); j++) {
            //if(_pdo[i][j]>0 && _pdo[i][j]<1024 && _chip[i]<=7) {
            if(_pdo[i][j]>0 && _pdo[i][j]<1024 && _threshold[i][j]==1) {
                if(_bcid[i][j]>0) {
                    has_hits = true;


                    ////////////////////////////////////////////
                    // RAW HISTOS
                    ////////////////////////////////////////////

                    // fill the raw charge histo
                    if(_chip[i] >=0 && _chip[i] < 8) {
                        h_pdo_raw.at(_chip[i])->Fill(_pdo[i][j]);
                    }

                    // fill vmm channel hit maps
                    if(_chip[i] >= 0 && _chip[i] < 8) {
                        h_vmm_hits.at(_chip[i])->Fill(_channel[i][j]);
                    }

                    int strip = TChamberMapping(_chip[i], _channel[i][j]);
                    if(_chip[i]>=0 && _chip[i]<4) {
                        h_chamber_hits.at(0)->Fill(strip);
                        h2_event_hit_map.at(0)->Fill(strip, _pdo[i][j]);
                    }
                    else if(_chip[i]>=4 && _chip[i]<8) {
                        h_chamber_hits.at(1)->Fill(strip);
                        h2_event_hit_map.at(1)->Fill(strip, _pdo[i][j]);
                    }

                    /////////////////////////////////////////////
                    // CLUSTER
                    /////////////////////////////////////////////
                    tmp_cluster_channel.push_back(_channel[i][j]);
                    tmp_cluster_strips.push_back(strip);
                    tmp_cluster_pdos.push_back(_pdo[i][j]);
                    tmp_cluster_bcids.push_back(_bcid[i][j]);
                    tmp_cluster_grays.push_back(_gray[i][j]);
                    tmp_cluster_chips.push_back(_chip[i]);

                } // bcid > 0
            } // within valid ranges
        } // j
    } // i

//    if(has_hits) {
//        ////////////////////////////////////////////////
//        // 2D strip vs. pdo hit map
//        ////////////////////////////////////////////////
//        for(int i = 0; i < 2; i++) {
//            c_event_hits->cd(i+1);
//            h2_event_hit_map.at(i)->Draw();
//            c_event_hits->Update();
//        }
//    }
//    getchar();
    
}

void VMMAna::sort_by_bcid()
{
    if(!(tmp_cluster_grays.size())) return;
    for(int i = 0; i < (int)tmp_cluster_grays.size(); i++) {
        for(int j = (int)tmp_cluster_grays.size()-1; j > i; j--) {
            if(tmp_cluster_grays.at(j-1) > tmp_cluster_grays.at(j))
            {
                swap(tmp_cluster_channel.at(j-1), tmp_cluster_channel.at(j));
                swap(tmp_cluster_strips.at(j-1),  tmp_cluster_strips.at(j));
                swap(tmp_cluster_pdos.at(j-1),    tmp_cluster_pdos.at(j));
                swap(tmp_cluster_bcids.at(j-1),   tmp_cluster_bcids.at(j));
                swap(tmp_cluster_chips.at(j-1),   tmp_cluster_chips.at(j));
                swap(tmp_cluster_grays.at(j-1),   tmp_cluster_grays.at(j));
            }
        } // j
    } // i
}

void VMMAna::performDuplicateRemoval()
{
    duplicate_found = false;
    n_duplicate_per_channel = 0;
    n_duplicate_per_event = 0;
    n_crossfire_per_event = 0;

    bool remove_first_duplicate = false; // ?? 
    if(tmp_cluster_strips.size()>1) {
        for(int ichan = 0; ichan < (int)tmp_cluster_strips.size(); ) {
            bool flip_i = false;
            n_duplicate_per_channel = 0;

            for(int jchan=ichan+1; jchan < (int)tmp_cluster_strips.size(); jchan++) {
                if( (tmp_cluster_strips.at(ichan) == tmp_cluster_strips.at(jchan)) &&
                    (tmp_cluster_pdos.at(ichan) == tmp_cluster_pdos.at(jchan)) &&
                    (tmp_cluster_bcids.at(ichan) == tmp_cluster_bcids.at(jchan)) ) {

                    // fill the duplicate histogram
                    h_duplicate_chan.at(tmp_cluster_chips.at(jchan))->Fill(tmp_cluster_channel.at(jchan));

                    if(tmp_cluster_chips.at(jchan) >=0 && tmp_cluster_chips.at(jchan)<4)
                        h_duplicate_strip.at(0)->Fill(tmp_cluster_strips.at(jchan));
                    else {
                        h_duplicate_strip.at(1)->Fill(tmp_cluster_strips.at(jchan));
                    }

                    // remove duplicate from cluster primitive
                    tmp_cluster_pdos.erase(tmp_cluster_pdos.begin() +       jchan);
                    tmp_cluster_channel.erase(tmp_cluster_channel.begin() + jchan);
                    tmp_cluster_strips.erase(tmp_cluster_strips.begin() +   jchan);
                    tmp_cluster_bcids.erase(tmp_cluster_bcids.begin() +     jchan);
                    tmp_cluster_grays.erase(tmp_cluster_grays.begin() +     jchan);
                    tmp_cluster_chips.erase(tmp_cluster_chips.begin() +     jchan);
                    jchan--;
                    duplicate_found = true;
                    n_duplicate_per_channel++;
                    n_duplicate_per_event++;

                } // duplicate found
                else if( (tmp_cluster_chips.at(ichan) == tmp_cluster_chips.at(jchan)) &&
                         (fabs(tmp_cluster_channel.at(ichan)-tmp_cluster_channel.at(jchan))==1) &&
                         (tmp_cluster_pdos.at(ichan) == tmp_cluster_pdos.at(jchan)) &&
                         (tmp_cluster_bcids.at(ichan) == tmp_cluster_bcids.at(jchan))) {

                        if(tmp_cluster_channel.at(ichan) > tmp_cluster_channel.at(jchan)) {
                            n_crossfire_per_event++;

                            tmp_cluster_pdos.erase(tmp_cluster_pdos.begin() + jchan);
                            tmp_cluster_channel.erase(tmp_cluster_channel.begin() + jchan);
                            tmp_cluster_strips.erase(tmp_cluster_strips.begin() + jchan);
                            tmp_cluster_bcids.erase(tmp_cluster_bcids.begin() + jchan);
                            tmp_cluster_grays.erase(tmp_cluster_grays.begin() + jchan);
                            tmp_cluster_chips.erase(tmp_cluster_chips.begin() + jchan);
                            jchan--;
                        } // i > j
                        else {
                            n_crossfire_per_event++;

                            tmp_cluster_pdos.erase(tmp_cluster_pdos.begin() +       ichan);
                            tmp_cluster_channel.erase(tmp_cluster_channel.begin() + ichan);
                            tmp_cluster_strips.erase(tmp_cluster_strips.begin() +   ichan);
                            tmp_cluster_bcids.erase(tmp_cluster_bcids.begin() +     ichan);
                            tmp_cluster_grays.erase(tmp_cluster_grays.begin() +     ichan);
                            tmp_cluster_chips.erase(tmp_cluster_chips.begin() +     ichan);
                            //ichan--;
                            flip_i = true;
                        }
                } // crossfire 
            } // jchan

            if(flip_i) ichan = ichan;
            else { ichan += 1; }

            // fill with total number of duplicates for this channel (ichan)
            h_total_duplicates->Fill(n_duplicate_per_channel);
        } // ichan
    } // >1 strip

    h_total_cross_duplicates->Fill(n_crossfire_per_event);
    h_total_duplicates->Fill(n_duplicate_per_event);

}

void VMMAna::draw_raw_histos()
{
    ////////////////////////////////////////////////
    // charge histos
    ////////////////////////////////////////////////
    c_charge = new TCanvas("c_charge", "", 1400, 800);

    c_charge->Divide(4,2);
    for(int i = 0; i < 8; i++) {
        c_charge->cd(i+1);
        h_pdo_raw.at(i)->Draw();
    }
    //c_charge->cd(2);
    //h_pdo_raw->Draw();

    ////////////////////////////////////////////////
    // hit maps (vmm)
    ////////////////////////////////////////////////
    c_hit_maps = new TCanvas("c_hit_map", "", 1000, 800);

    c_hit_maps->Divide(4,2);
    for(int i = 0; i < 8; i++) {
        c_hit_maps->cd(i+1);
        h_vmm_hits.at(i)->Draw();
    }

    ////////////////////////////////////////////////
    // hit maps (T-chamber)
    ////////////////////////////////////////////////
    c_chamber_hits = new TCanvas("c_chamber_hits", "", 1000, 800);
    c_chamber_hits->Divide(1,2);
    for(int i = 0; i < 2; i++) {
        c_chamber_hits->cd(i+1);
        h_chamber_hits.at(i)->Draw();
    }
}

void VMMAna::draw_duplicate_histos()
{
    /////////////////////////////////////////
    // duplicates per VMM channel (8)
    /////////////////////////////////////////
    c_duplicate_channel = new TCanvas("c_duplicate_channel", "", 1000, 800);
    c_duplicate_channel->Divide(4,2);
    for(int i = 0; i < 8; i++) {
        c_duplicate_channel->cd(i+1);
        h_duplicate_chan.at(i)->Draw();
    }

    /////////////////////////////////////////
    // duplicates per T-chamber strip (2)
    /////////////////////////////////////////
    c_duplicate_strip = new TCanvas("c_duplicate_strip", "", 1000, 800);
    c_duplicate_strip->Divide(1,2);
    for(int i = 0; i < 2; i++) {
        c_duplicate_strip->cd(i+1);
        h_duplicate_strip.at(i)->Draw();
    }

    /////////////////////////////////////////
    // total duplicates
    /////////////////////////////////////////
    c_total_duplicates = new TCanvas("c_total_duplicates", "", 1000, 500);
    c_total_duplicates->Divide(2,1);
    c_total_duplicates->cd(1);
    h_total_duplicates->Draw();
    c_total_duplicates->cd(2);
    h_total_cross_duplicates->Draw();

}

int VMMAna::TChamberMapping(int vmm_id, int vmm_channel)
{
    int strip = 0;
    if(vmm_id%2==0) {
        if(vmm_channel== 0 )   strip= 32 ;
        if(vmm_channel== 1 )   strip= 33 ;
        if(vmm_channel== 2 )   strip= 31 ;
        if(vmm_channel== 3 )   strip= 34 ;
        if(vmm_channel== 4 )   strip= 30 ;
        if(vmm_channel== 5 )   strip= 35 ;
        if(vmm_channel== 6 )   strip= 29 ;
        if(vmm_channel== 7 )   strip= 36 ;
        if(vmm_channel== 8 )   strip= 28 ;
        if(vmm_channel== 9 )   strip= 37 ;
        if(vmm_channel== 10)   strip= 27 ;
        if(vmm_channel== 11)   strip= 38 ;
        if(vmm_channel== 12)   strip= 26 ;
        if(vmm_channel== 13)   strip= 39 ;
        if(vmm_channel== 14)   strip= 25 ;
        if(vmm_channel== 15)   strip= 40 ;
        if(vmm_channel== 16)   strip= 24 ;
        if(vmm_channel== 17)   strip= 41 ;
        if(vmm_channel== 18)   strip= 23 ;
        if(vmm_channel== 19)   strip= 42 ;
        if(vmm_channel== 20)   strip= 22 ;
        if(vmm_channel== 21)   strip= 43 ;
        if(vmm_channel== 22)   strip= 21 ;
        if(vmm_channel== 23)   strip= 44 ;
        if(vmm_channel== 24)   strip= 20 ;
        if(vmm_channel== 25)   strip= 45 ;
        if(vmm_channel== 26)   strip= 19 ;
        if(vmm_channel== 27)   strip= 46 ;
        if(vmm_channel== 28)   strip= 18 ;
        if(vmm_channel== 29)   strip= 47 ;
        if(vmm_channel== 30)   strip= 17 ;
        if(vmm_channel== 31)   strip= 48 ;
        if(vmm_channel== 32)   strip= 16 ;
        if(vmm_channel== 33)   strip= 49 ;
        if(vmm_channel== 34)   strip= 15 ;
        if(vmm_channel== 35)   strip= 50 ;
        if(vmm_channel== 36)   strip= 14 ;
        if(vmm_channel== 37)   strip= 51 ;
        if(vmm_channel== 38)   strip= 13 ;
        if(vmm_channel== 39)   strip= 52 ;
        if(vmm_channel== 40)   strip= 12 ;
        if(vmm_channel== 41)   strip= 53 ;
        if(vmm_channel== 42)   strip= 11 ;
        if(vmm_channel== 43)   strip= 54 ;
        if(vmm_channel== 44)   strip= 10 ;
        if(vmm_channel== 45)   strip= 55 ;
        if(vmm_channel== 46)   strip= 9  ;
        if(vmm_channel== 47)   strip= 56 ;
        if(vmm_channel== 48)   strip= 8  ;
        if(vmm_channel== 49)   strip= 57 ;
        if(vmm_channel== 50)   strip= 7  ;
        if(vmm_channel== 51)   strip= 58 ;
        if(vmm_channel== 52)   strip= 6  ;
        if(vmm_channel== 53)   strip= 59 ;
        if(vmm_channel== 54)   strip= 5  ;
        if(vmm_channel== 55)   strip= 60 ;
        if(vmm_channel== 56)   strip= 4  ;
        if(vmm_channel== 57)   strip= 61 ;
        if(vmm_channel== 58)   strip= 3  ;
        if(vmm_channel== 59)   strip= 62 ;
        if(vmm_channel== 60)   strip= 2  ;
        if(vmm_channel== 61)   strip= 63 ;
        if(vmm_channel== 62)   strip= 1  ;
        if(vmm_channel== 63)   strip= 64 ;
    }
    else {
        if(vmm_channel== 0 )   strip=   128  ;
        if(vmm_channel== 1 )   strip=   65   ;
        if(vmm_channel== 2 )   strip=   127  ;
        if(vmm_channel== 3 )   strip=   66   ;
        if(vmm_channel== 4 )   strip=   126  ;
        if(vmm_channel== 5 )   strip=   67   ;
        if(vmm_channel== 6 )   strip=   125  ;
        if(vmm_channel== 7 )   strip=   68   ;
        if(vmm_channel== 8 )   strip=   124  ;
        if(vmm_channel== 9 )   strip=   69   ;
        if(vmm_channel== 10)   strip=   123  ;
        if(vmm_channel== 11)   strip=   70   ;
        if(vmm_channel== 12)   strip=   122  ;
        if(vmm_channel== 13)   strip=   71   ;
        if(vmm_channel== 14)   strip=   121  ;
        if(vmm_channel== 15)   strip=   72   ;
        if(vmm_channel== 16)   strip=   120  ;
        if(vmm_channel== 17)   strip=   73   ;
        if(vmm_channel== 18)   strip=   119  ;
        if(vmm_channel== 19)   strip=   74   ;
        if(vmm_channel== 20)   strip=   118  ;
        if(vmm_channel== 21)   strip=   75   ;
        if(vmm_channel== 22)   strip=   117  ;
        if(vmm_channel== 23)   strip=   76   ;
        if(vmm_channel== 24)   strip=   116  ;
        if(vmm_channel== 25)   strip=   77   ;
        if(vmm_channel== 26)   strip=   115  ;
        if(vmm_channel== 27)   strip=   78   ;
        if(vmm_channel== 28)   strip=   114  ;
        if(vmm_channel== 29)   strip=   79   ;
        if(vmm_channel== 30)   strip=   113  ;
        if(vmm_channel== 31)   strip=   80   ;
        if(vmm_channel== 32)   strip=   112  ;
        if(vmm_channel== 33)   strip=   81   ;
        if(vmm_channel== 34)   strip=   111  ;
        if(vmm_channel== 35)   strip=   82   ;
        if(vmm_channel== 36)   strip=   110  ;
        if(vmm_channel== 37)   strip=   83   ;
        if(vmm_channel== 38)   strip=   109  ;
        if(vmm_channel== 39)   strip=   84   ;
        if(vmm_channel== 40)   strip=   108  ;
        if(vmm_channel== 41)   strip=   85   ;
        if(vmm_channel== 42)   strip=   107  ;
        if(vmm_channel== 43)   strip=   86   ;
        if(vmm_channel== 44)   strip=   106  ;
        if(vmm_channel== 45)   strip=   87   ;
        if(vmm_channel== 46)   strip=   105  ;
        if(vmm_channel== 47)   strip=   88   ;
        if(vmm_channel== 48)   strip=   104  ;
        if(vmm_channel== 49)   strip=   89   ;
        if(vmm_channel== 50)   strip=   103  ;
        if(vmm_channel== 51)   strip=   90   ;
        if(vmm_channel== 52)   strip=   102  ;
        if(vmm_channel== 53)   strip=   91   ;
        if(vmm_channel== 54)   strip=   101  ;
        if(vmm_channel== 55)   strip=   92   ;
        if(vmm_channel== 56)   strip=   100  ;
        if(vmm_channel== 57)   strip=   93   ;
        if(vmm_channel== 58)   strip=   9    ;
        if(vmm_channel== 59)   strip=   94   ;
        if(vmm_channel== 60)   strip=   8    ;
        if(vmm_channel== 61)   strip=   95   ;
        if(vmm_channel== 62)   strip=   7    ;
        if(vmm_channel== 63)   strip=   96   ;
    }

    if(vmm_id%2==0 && (vmm_id==2 || vmm_id==3)) {
        strip+=128;
    }
    else if(!(vmm_id%2==0) && (vmm_id==6 || vmm_id==7)) {
        strip+=128;
    }

    return strip;

}
