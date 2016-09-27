#include "t6t7_ana.hh"

//std/stl
#include <iostream>
#include <string>
#include <set>

//ROOT
#include "TSystem.h"
#include "TRandom.h"
#include "TF1.h"

using namespace t6t7;
using namespace std;

double PITCH = 0.4;
double GAIN = 3.0;
int MIN_CLUSTER_SIZE = 3;

//--------------------------------------------------------------//
// Constructor
VMMAna::VMMAna() :
    n_total(0),
    n_T6offset_counter(0),
    m_chain(0),
    m_dbg(false),
    m_calib_file(""),
    m_calib(false)
{
    gSystem->Load("/Users/dantrim/workarea/NSW/myanalysis/vmm-analysis/vector_lib/libMylib");
}
//--------------------------------------------------------------//
// comparison operators for hits/clusters
struct bcidSmaller { // will sort to have earlier hits first (use the gray decoded BCID!)
    bool operator()(vmm::Hit a, vmm::Hit b) { return a.gray() < b.gray(); }
} byBCID;

struct clusPdoLarger {
    bool operator()(vmm::Cluster a, vmm::Cluster b) { return a.pdo() > b.pdo(); }
} byClusPdo;

struct stripSmaller {
    bool operator()(vmm::Hit a, vmm::Hit b) { return a.strip() < b.strip(); }
} byStrip;

struct clusPdoEqual {
    bool operator()(vmm::Cluster a, vmm::Cluster b) { return ((a.pdo() == b.pdo()) && (a.chamber()==b.chamber())); }
} clusPdoEqual;


//--------------------------------------------------------------//
// TSelector::Begin
void VMMAna::Begin(TTree* /*tree*/)
{
    
    cout << " ------------------------ VMMAna::Begin ----------------------- " << endl;

}
//--------------------------------------------------------------//
// TSelector::SlaveBegin
void VMMAna::SlaveBegin(TTree* /*tree*/)
{
    return;
}
//--------------------------------------------------------------//
// TSelector::Init
void VMMAna::Init(TTree* /*tree*/)
{
    initializeNtupleVariablesAndBranches();

    initializeHistograms();
}
//--------------------------------------------------------------//
// TSelector::Notify
Bool_t VMMAna::Notify()
{
    return kTRUE;
}

//--------------------------------------------------------------//
void VMMAna::loadCalibration()
{
    if(calib()) {
        cout << "VMMAna::loadCalibration    Handling of calibration data is not supported!" << endl;
        cout << "VMMAna::loadCalibration     >>> Exitting." << endl;
        exit(1);
    }
    else {
        // T6
        std::vector<double> t6_channel;
        std::vector<double> t6_gain;
        std::vector<double> t6_pedestal;
        for(int i = 0; i < 256; i++) {
            t6_channel.push_back(i);
            t6_gain.push_back(3.);
            t6_pedestal.push_back(32.);
        } // i
        m_channel_calib.push_back(t6_channel);
        m_channel_gain_calib.push_back(t6_gain);
        m_channel_pedestal_calib.push_back(t6_pedestal);
        // T7
        std::vector<double> t7_channel;
        std::vector<double> t7_gain;
        std::vector<double> t7_pedestal;
        for(int i = 0; i < 256; i++) {
            t7_channel.push_back(i);
            t7_gain.push_back(3.);
            t7_pedestal.push_back(32.);
        } // i
        m_channel_calib.push_back(t7_channel);
        m_channel_gain_calib.push_back(t7_gain);
        m_channel_pedestal_calib.push_back(t7_pedestal);
    }

}
//--------------------------------------------------------------//
// TSelector::Process
Bool_t VMMAna::Process(Long64_t entry)
{
    static Long64_t chainEntry = -1;
    chainEntry++;
    m_chain->GetEntry(chainEntry);

    if(chainEntry%10000==0 || m_dbg) {
        cout << " ***** Processing entry " << chainEntry << endl;
    }

    loadCalibration();


    // offset in T6 strips by +2 w.r.t. T7 is done here
    int n_hits = loadChamberHits();

    fillRawHistograms();

    // sort the hits according to larger bcid (last will be latest)
    sort(hits.begin(), hits.end(), byBCID);

    #warning not performing duplicate removal

    removeDuplicateHits();

    //cout << "!! TOTAL OF " << n_hits << " HITS !! " << endl;
    //cout << "print hits: " << endl;
    //for(int i = 0; i < (int)hits.size(); i++) {
    //    hits.at(i).print();
    //}


    // make preliminary clusters
    int n_clusters = makeClusters();
    getClusterChargeAndPosition();

    removeDuplicateClusters();
    removeDuplicateStripsInClusters();
    int n_clusters_after_overlap = clusters.size();

    fillChargeHistos();

    fillClusterMultiplicityHistos();

    fillClusterHistos();

    fillLeadStripChargeHistos();

    // require certain number of strips per cluster
    if(!(passClusterSize(MIN_CLUSTER_SIZE))) return false;


    // fill Q ratio histos
    //fillChargeRatioHistos();

    // now apply some cuts and apply matching between T6/T7
    // > here require exactly 1 cluster per layer
    if(!(nClusterPerLayer(1, true))) return false;

    // check if clusters are overlapping, require X number of
    // strips to be overlapping between T6 and T7
    if(!(overlappingClusters(1))) return false;

    fillCorrHistos();

    n_total++;

    fillClusterHistos_OR();
    

    return true;
}
//--------------------------------------------------------------//
void VMMAna::fillRawHistograms()
{
    for(auto& hit : hits) {
        if(hit.chip()>=8 || hit.chip()<0) continue;
        if(!(hit.bcid()>0)) continue;
        if(hit.pdo()>0 && hit.pdo()<1024 && hit.threshold()==1) {

            int chip_no = hit.chip();

            // raw pdo histos
            h_pdo_raw.at(chip_no)->Fill(hit.pdo());

            // raw vmm channel hits
            h_vmm_hits_raw.at(chip_no)->Fill(hit.channel());

            // raw T-chamber hit map
            if(chip_no>=0 && chip_no<4) {
                h_strip_hits_raw.at(0)->Fill(hit.strip());
            }
            else {
                h_strip_hits_raw.at(1)->Fill(hit.strip());
            }

        }
    } // hit

}

//--------------------------------------------------------------//
void VMMAna::fillChargeHistos()
{
    if(hits.size()>=1) {
        for(auto hit : hits) {
            int pdo = hit.pdo();
            double pdo_calib = 0.0;

            if(calib()) {
                cout << "VMMAna::fillChargeHistos    Attempting to apply calibration. This is "
                    << "not setup. Exitting." << endl;
                exit(1);
            }
            else {
                pdo_calib = (GAIN*pdo) / m_channel_gain_calib[hit.chamber()][hit.channel()];
            }
            if(pdo>0) {
                h_pdo.at(hit.chamber())->Fill(pdo);
                h_pdo_calib.at(hit.chamber())->Fill(pdo_calib);
            }
        }
    }
    
}
//--------------------------------------------------------------//
void VMMAna::fillClusterMultiplicityHistos()
{
    int n_cluster_t6 = 0;
    int n_cluster_t7 = 0;
    for(auto cluster : clusters) {
        if(cluster.chamber()==0) n_cluster_t6++;
        else if(cluster.chamber()==1) n_cluster_t7++;
    }
    h_num_cluster.at(0)->Fill(n_cluster_t6);
    h_num_cluster.at(1)->Fill(n_cluster_t7);

    if(clusters.size()>0) {
        for(auto& cluster : clusters) {
            if(cluster.chamber()==0) {
                h_cluster_size.at(0)->Fill(cluster.size());
                h2_num_cluster_vs_cluster_size.at(0)->Fill(n_cluster_t6, cluster.size());
            } // T6
            else if(cluster.chamber()==1) {
                h_cluster_size.at(1)->Fill(cluster.size());
                h2_num_cluster_vs_cluster_size.at(1)->Fill(n_cluster_t7, cluster.size());
            }
        }
    }
}
//--------------------------------------------------------------//
std::string VMMAna::intToBinaryStr(int number)
{
    string bin;
    char holder = ' ';
    while(number!=0) {
        holder=number%2+'0';
        bin=holder+bin;
        number/=2;
    }
    return bin;
}
int VMMAna::binaryStrToInt(string binary)
{
    char binaryChar[binary.length()];
    strcpy(binaryChar, binary.c_str());
    char *ptr;
    int parsed = strtol(binaryChar, &ptr, 2);
    return parsed;
}
//--------------------------------------------------------------//
void VMMAna::fillCorrHistos()
{
    // correlation histos
    if(clusters.size()>0) {
        double t6_cluster_charge = 0.0;
        double t7_cluster_charge = 0.0;

        vector<int> t6_strips;
        vector<int> t7_strips;

        for(auto& cluster : clusters) {
            if(cluster.size()>=1) {

             
                int cluster_size_cut = 2;
                if(!(cluster.pdo()>0.)) continue;
                if(!(cluster.size()>=cluster_size_cut)) continue;

                if(cluster.chamber()==0) {
                    t6_cluster_charge += cluster.pdoCalibrated();

                    for(auto& strip : cluster.hits()) {
                        t6_strips.push_back(strip.strip());
                    }

                }
                else if(cluster.chamber()==1) {
                    t7_cluster_charge += cluster.pdoCalibrated();

                    for(auto& strip : cluster.hits()) {
                        t7_strips.push_back(strip.strip());
                    }
                }

            }
        } // cluster
        if(t6_cluster_charge>0 || t7_cluster_charge>0)
            h2_cl_charge_T6vsT7->Fill(t6_cluster_charge, t7_cluster_charge);

        if(t6_strips.size() != t7_strips.size()) {
            if(t6_strips.size() < t7_strips.size()) {
                int n = (int)(t7_strips.size() - t6_strips.size());
                for(int i = 0; i < n; i++) t6_strips.push_back(-1);
            }
            else if(t7_strips.size() < t6_strips.size()) {
                int n = (int)(t6_strips.size() - t7_strips.size());
                for(int i = 0; i < n; i++) t7_strips.push_back(-1);
            }
        }

        if(t6_strips.size()!=t7_strips.size()) {
            cout << "BLAH" << endl;
            cout << "   t6: " << t6_strips.size() <<endl;
            cout << "   t7: " << t7_strips.size() << endl;
            exit(1);
        }

        for(int i = 0; i < (int)t6_strips.size(); i++) {
            h2_cl_strip_T6vsT7->Fill(t6_strips.at(i), t7_strips.at(i));
        }
    }

}
//--------------------------------------------------------------//
void VMMAna::fillClusterHistos()
{

    if(clusters.size()>0) {
        for(auto& cluster : clusters) {
            if(cluster.size()>=1) {

                int cluster_size_cut = 2; // don't fill unless clusters are >=2 strips
                if(!(cluster.pdo()>0.)) continue;
                if(!(cluster.size()>=cluster_size_cut)) continue;

                // charge histos for cluster
                h_cl_charge.at(cluster.chamber())->Fill(cluster.pdo());
                h_cl_charge_calib.at(cluster.chamber())->Fill(cluster.pdoCalibrated());

                // position histos for cluster
                h_cl_position.at(cluster.chamber())->Fill(cluster.positionCentroid());
                // currently not doing calibration, so don't fill calibrated position
            } // cl size >=1
        } // cluster
    }
}
//--------------------------------------------------------------//
void VMMAna::fillClusterHistos_OR()
{
    if(clusters.size()>0) {
        for(auto& cluster : clusters) {
            if(cluster.size()>=1) {
                int cluster_size_cut = 2;
                if(!(cluster.pdo()>0.)) continue;
                if(!(cluster.size()>=cluster_size_cut)) continue;

                // charge histos
                h_cl_charge_OR.at(cluster.chamber())->Fill(cluster.pdo());
                h_cl_charge_calib_OR.at(cluster.chamber())->Fill(cluster.pdoCalibrated());

                // position histos
                h_cl_position_OR.at(cluster.chamber())->Fill(cluster.positionCentroid());

            } // cluster size >= 1
        } // cluster
    } // size > 0
}
//--------------------------------------------------------------//
void VMMAna::fillChargeRatioHistos()
{
//    cout << "h ratio size : " << h_cl_strip_qratio.size() << endl;
    for(auto cluster : clusters) {
        vector<vmm::Hit> strips = cluster.hits();

        vector<vmm::Hit> strips_sorted_by_strip;
        for(auto st : strips) strips_sorted_by_strip.push_back(st);
        std::sort(strips_sorted_by_strip.begin(), strips_sorted_by_strip.end(), byStrip);

        int n_strips = (int)strips_sorted_by_strip.size();
        //for(auto& st : strips_sorted_by_strip) {
        //    cout << st.strip() << "  " << st.pdo() << endl;
        //}
        //cout << "-------------------" << endl;
        for(int i = 0; i < n_strips; i++) {
            if(i==(n_strips-1)) continue;
            if(strips_sorted_by_strip.at(i+1).pdo()>0) {
                double ratio = strips_sorted_by_strip.at(i).pdo() / strips_sorted_by_strip.at(i+1).pdo() * 1.0;
                h_cl_strip_qratio.at(i)->Fill(ratio);
            }
        } // i
    }
}

//--------------------------------------------------------------//
void VMMAna::fillLeadStripChargeHistos()
{
    for(auto cluster : clusters) {
        int pdo_largest = -1;
        for(auto strip : cluster.hits()) {
            if(strip.pdo() > pdo_largest) pdo_largest = strip.pdo();
        } // strip
        if(pdo_largest>0)
            h_cl_lead_stripq.at(cluster.chamber())->Fill(pdo_largest);
    } // cluster
}
//--------------------------------------------------------------//
bool VMMAna::passClusterSize(int n_strip_per_cluster)
{
    vector<vmm::Cluster> tmpClusters;

    for(auto cluster : clusters) {
        if(cluster.size() >= n_strip_per_cluster) {
            tmpClusters.push_back(cluster);
        }
    }
    clusters.clear();

    for(auto cluster : tmpClusters) {
        clusters.push_back(cluster);
    }
    return clusters.size();

}
//--------------------------------------------------------------//
bool VMMAna::nClusterPerLayer(int n_per_layer, bool exclusive)
{
    bool pass = false;

    int n_t6 = 0;
    int n_t7 = 0;

    for(auto cluster : clusters) {
        if(cluster.chamber()==0) n_t6++;
        else if(cluster.chamber()==1) n_t7++;
    } // cluster

    if(exclusive) {
        if(n_t6==n_per_layer && n_t7==n_per_layer) pass = true;
    }
    else if(!exclusive) {
        if(n_t6>=n_per_layer && n_t7>=n_per_layer) pass = true;
    }
    return pass;

}
//--------------------------------------------------------------//
bool VMMAna::overlappingClusters(int n_min_common)
{

    bool overlap = false;

    // actually lets do overlap based on centroid position
    double t6_position = 0.0;
    double t7_position = 0.0;

    if(clusters.size()==2) {
        for(auto& cluster : clusters) {
            // T6
            if(cluster.chamber()==0) {
                t6_position = cluster.positionCentroidCalibrated();
            }
            // T7
            else if(cluster.chamber()==1) {
                t7_position = cluster.positionCentroidCalibrated();
            }
        }
    }
    else {
        cout << "VMMAna::overlappingClusters    ERROR Can only handle case with 1 cluster per "
            << "layer and you have more than 2 clusters! Exitting." << endl;
        exit(1);
    }

    //h_cl_residual->Fill(residual);

    // find the offset of T7
    // correct T7 position 
    //t7_position = 1.06 * t7_position;
    //t7_position *= t7_position/t6_position;
    double residual = t7_position - t6_position;
    h2_cl_position_T6vsT7->Fill(t6_position, t7_position);

    //double overlap_condition = MIN_CLUSTER_SIZE*PITCH + PITCH;
    //double overlap_condition = MIN_CLUSTER_SIZE*1.25;
    double overlap_condition = 0.1 * 5; // 5 sigma resolution, assuming 100 um resolution (0.1 is in mm)
    //cout << "t6_position - t7_position > " << t6_position << " - " << t7_position << " = " << (t6_position-t7_position) << "  residual: " << residual << "   OR condition: " << overlap_condition << endl;
    if(fabs(t7_position - t6_position) <= overlap_condition) overlap = true;

    h_cl_residual->Fill(residual);
    //if(overlap) h_cl_residual->Fill(residual);

    return overlap;

     
/*
    // gather the strip numbers associated with each of the clusters in T6/T7

    vector<int> strips_t6;
    vector<int> strips_t7;

    strips_t6.clear();
    strips_t7.clear(); 

    for(auto cluster : clusters) {
        if(cluster.chamber()==0) {
            for(auto strip : cluster.hits()) {
                strips_t6.push_back(strip.strip());
            }
        }
        else if(cluster.chamber()==1) {
            for(auto strip : cluster.hits()) {
                strips_t7.push_back(strip.strip());
            }
        }
    }

    //for(auto strip_no : strips_t6) {
    //    cout << strip_no << "  ";
    //}
    //cout << endl;
    //cout << "----------------------------------" << endl;

    // now check that T6 strips are overlapping with T7
    // > will require >=n_min_common strips in common between the chambers
    int n_overlap = 0;
    for(int i6 = 0; i6 < (int)strips_t6.size(); i6++) {
        int t6Strip = strips_t6.at(i6);
        for(int i7 = 0; i7 < (int)strips_t7.size(); i7++) {
            int t7Strip = strips_t7.at(i7);
            if(t6Strip==t7Strip) {
                n_overlap++;
                break; // move to next strip in T6
            }
        }
    } // loop over T6 strip numbers

    if(n_overlap >= n_min_common) overlap = true;
    else { overlap = false; }
    if(overlap) h_cl_residual->Fill(residual);
    
    return overlap;
*/
}
//--------------------------------------------------------------//
// initialize histograms
void VMMAna::initializeHistograms()
{
    initRawHistos();

    initClusterMultiplicityHistos();

    initCleanedChargeHistos();

    initClusterChargeHistos();

    initClusterPositionHistos();

    initResidualHistos();

    initCorrHistos();

    initChargeRatioHistos();

    initLeadStripChargeHistos();

}
void VMMAna::initRawHistos()
{

    // raw pdo straight from VMM over all channels
    stringstream name;
    for(int i = 0; i < 8; i++) {
        name.str("");
        name << "VMM " << i << " Charge Distribution (all chan, no calib)";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 400, 0, 1200);
        h->GetXaxis()->SetTitle("PDO [ADC counts]");
        h->GetYaxis()->SetTitle("# of Entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_pdo_raw.push_back(h); 
    } // i
    name.str("");

    // raw vmm channel hits (no cleaning)
    for(int i = 0; i < 8; i++) {
        name << "h_vmm_chan_hits_" << i;
        TH1F* h = new TH1F(name.str().c_str(), "", 63,0,63); name.str("");
        name << "VMM " << i << " Raw Channel Hits";
        h->SetTitle(name.str().c_str()); name.str("");
        name << "VMM " << i << " channel";
        h->GetXaxis()->SetTitle(name.str().c_str()); name.str("");
        h->GetYaxis()->SetTitle("# of Entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_vmm_hits_raw.push_back(h);
    } // i
    name.str("");

    // raw chamber hits (no cleaning)
    for(int i = 0; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " Chamber Raw Hit Map";
        if(i==1) name << " (re-mapped to T6)";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 256, 0, 256);
        h->GetXaxis()->SetTitle("MicroMegas Strip");
        h->GetYaxis()->SetTitle("# of Entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_strip_hits_raw.push_back(h);
    } // i
    name.str("");

}
void VMMAna::initClusterMultiplicityHistos()
{
    // number of clusters in each chamber
    stringstream name;
    for(int i = 0; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Number of Clusters/Event";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 20, 0, 20);
        h->GetXaxis()->SetTitle("# of clusters");
        h->GetYaxis()->SetTitle("# of entries");
        h_num_cluster.push_back(h);
    }
    name.str("");

    // cluster size at each chamber
    for(int i = 0; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Size of Clusters";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 10, 0, 10);
        h->GetXaxis()->SetTitle("# of strips");
        h->GetYaxis()->SetTitle("# of Entries");
        h_cluster_size.push_back(h);
    }
    name.str("");

    // cl size vs multiplicity
    for(int i = 0; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Cluster Mult. vs. Size";
        TH2F* h = new TH2F(name.str().c_str(), name.str().c_str(), 20, 0, 20, 10, 0, 10);
        h->GetXaxis()->SetTitle("# of clusters per event");
        h->GetYaxis()->SetTitle("# of strips per cluster");
        h2_num_cluster_vs_cluster_size.push_back(h);
    }
    name.str("");


    //h_num_cluster = new TH1F("h_num_cluster", "Number of Clusters Per Event", 20, 0, 20);

    //h_cluster_size = new TH1F("h_cluster_size", "Size of Clusters", 10, 0, 10);

    //h2_num_cluster_vs_cluster_size = new TH2F("num_cluster_vs_cluster_size",
    //        "Cluster Multiplicity vs. Cluster Size", 20, 0, 20, 10, 0, 10);
    //h2_num_cluster_vs_cluster_size->GetXaxis()->SetTitle("# of clusters per event");
    //h2_num_cluster_vs_cluster_size->GetYaxis()->SetTitle("# of strips per cluster");

}
void VMMAna::initCleanedChargeHistos()
{
    stringstream name;

    // pdo of T-chamber strip
    for(int i = 0; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Charge distribution (all chan, no cleaning)";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 400, 0, 1200);
        h->GetXaxis()->SetTitle("PDO [ADC counts]");
        h->GetYaxis()->SetTitle("# of entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_pdo.push_back(h);
    }

 //   // pdo with ADC fix (smearing)
 //   for(int i = 0; i < 2; i++) {
 //       name.str("");
 //       name << "T" << (i==0 ? "6" : "7") << " : Charge distribution (all chan, w/ smearing)";
 //       TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 400, 0, 1200);
 //       h->GetXaxis()->SetTitle("PDO [ADC counts]");
 //       h->GetYaxis()->SetTitle("# of entries");
 //       h->GetYaxis()->SetTitleOffset(1.5);
 //       h_pdo_adc_fix.push_back(h);
 //   }

    // pdo with calibration
    for(int i = 0; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Charge distribution (all chan, w/ calibration)";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 400, 0, 1200);
        h->GetXaxis()->SetTitle("PDO [ADC counts]");
        h->GetYaxis()->SetTitle("# of entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_pdo_calib.push_back(h);
    }


}
void VMMAna::initClusterChargeHistos()
{
    stringstream name;

    // cluster charge
    for(int i = 0 ; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Cluster charge";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 200, 0, 4000);
        h->GetXaxis()->SetTitle("PDO [ADC counts]");
        h->GetYaxis()->SetTitle("# of entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_cl_charge.push_back(h);
    }

    // cluster charge (calib)
    for(int i = 0 ; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Cluster charge with Calibration";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 200, 0, 4000);
        h->GetXaxis()->SetTitle("PDO [ADC counts]");
        h->GetYaxis()->SetTitle("# of entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_cl_charge_calib.push_back(h);
    }

 //   // cluster charge (adc fix)
 //   for(int i = 0 ; i < 2; i++) {
 //       name.str("");
 //       name << "T" << (i==0 ? "6" : "7") << " : Cluster charge with Smearing";
 //       TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 200, 0, 4000);
 //       h->GetXaxis()->SetTitle("PDO [ADC counts]");
 //       h->GetYaxis()->SetTitle("# of entries");
 //       h->GetYaxis()->SetTitleOffset(1.5);
 //       h_cl_charge_adc_fix.push_back(h);
 //   }
    

    // --------------------------------------------------------------------- //
    //  POST OR
    // --------------------------------------------------------------------- //

    // cluster charge
    for(int i = 0 ; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Cluster charge (POST OR)";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 200, 0, 4000);
        h->GetXaxis()->SetTitle("PDO [ADC counts]");
        h->GetYaxis()->SetTitle("# of entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_cl_charge_OR.push_back(h);
    }
    
    // cluster charge (calib)
    for(int i = 0 ; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Cluster charge with Calibration (POST OR)";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 200, 0, 4000);
        h->GetXaxis()->SetTitle("PDO [ADC counts]");
        h->GetYaxis()->SetTitle("# of entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_cl_charge_calib_OR.push_back(h);
    }

 //   // cluster charge (adc fix)
 //   for(int i = 0 ; i < 2; i++) {
 //       name.str("");
 //       name << "T" << (i==0 ? "6" : "7") << " : Cluster charge with Smearing (POST OR)";
 //       TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 200, 0, 4000);
 //       h->GetXaxis()->SetTitle("PDO [ADC counts]");
 //       h->GetYaxis()->SetTitle("# of entries");
 //       h->GetYaxis()->SetTitleOffset(1.5);
 //       h_cl_charge_adc_fix_OR.push_back(h);
 //   }


}
void VMMAna::initClusterPositionHistos()
{
    stringstream name;
    // cluster position in T-chamber
    for(int i = 0; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Cluster Position in Chamber";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 200, 0, 256*0.4);
        h->GetXaxis()->SetTitle("Cluster Position [mm]");
        h->GetYaxis()->SetTitle("# of Entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_cl_position.push_back(h);
    }

    name.str("");

    for(int i = 0; i < 2; i++) {
        name.str("");
        name << "T" << (i==0 ? "6" : "7") << " : Cluster Position in Chamber (POST OR);";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 200, 0, 256*0.4);
        h->GetXaxis()->SetTitle("Cluster Position [mm]");
        h->GetYaxis()->SetTitle("# of Entries");
        h->GetYaxis()->SetTitleOffset(1.5);
        h_cl_position_OR.push_back(h);
    }
}
//--------------------------------------------------------------//
void VMMAna::initResidualHistos()
{
    stringstream name;
    name << "T6-T7 Cluster Hit Position Difference";
    h_cl_residual = new TH1F(name.str().c_str(), name.str().c_str(), 80, -2, 2);
    h_cl_residual->GetXaxis()->SetTitle("Residual [mm]");
    h_cl_residual->GetYaxis()->SetTitle("# of Entries");
    h_cl_residual->GetYaxis()->SetTitleOffset(1.5);

    name.str("");
    name << "T6 Cluster Position vs. T7 Cluster Position";
    h2_cl_position_T6vsT7 = new TH2F(name.str().c_str(), name.str().c_str(), 60, 0, 30, 60, 0, 30);
    h2_cl_position_T6vsT7->GetXaxis()->SetTitle("T6 Cluster Position [mm]");
    h2_cl_position_T6vsT7->GetYaxis()->SetTitle("T7 Cluster Position [mm]");

}
//--------------------------------------------------------------//
void VMMAna::initCorrHistos()
{
    stringstream name;
    name << "T6-T7 Cluster Charge Correlation";
    h2_cl_charge_T6vsT7 = new TH2F(name.str().c_str(), name.str().c_str(), 200, 0, 3000, 200, 0, 3000);
    h2_cl_charge_T6vsT7->GetXaxis()->SetTitle("T6 Cluster Charge [PDO counts]"); 
    h2_cl_charge_T6vsT7->GetYaxis()->SetTitle("T7 Cluster Charge [PDO counts]");


    name.str("");
    name << "T6-T7 Strip Hit Correlation";
    h2_cl_strip_T6vsT7 = new TH2F(name.str().c_str(), name.str().c_str(), 70, 0, 70, 70, 0, 70);
    h2_cl_strip_T6vsT7->GetXaxis()->SetTitle("T6 Strip Hits");
    h2_cl_strip_T6vsT7->GetYaxis()->SetTitle("T7 Strip Hits");

}
//--------------------------------------------------------------//
void VMMAna::initChargeRatioHistos()
{
    string num;
    string den;
    stringstream name;
    for(int i = 0; i < 10; i++) {
        num = i;
        den = i+1; 
        name << "Q_" << num << "/Q_" << den << " Ratio";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 20, 0, 1);
        h->GetXaxis()->SetTitle(name.str().c_str());
        h->GetYaxis()->SetTitle("# of Entries");
        h_cl_strip_qratio.push_back(h);
        name.str("");
    }
}
//--------------------------------------------------------------//
void VMMAna::initLeadStripChargeHistos()
{
    stringstream name;
    for(int i = 0; i < 2; i++) {
        name << "T" << (i==0 ? "6" : "7") << " Lead Strip Charge";
        TH1F* h = new TH1F(name.str().c_str(), name.str().c_str(), 200, 0, 1200);
        h->GetXaxis()->SetTitle("PDO [ADC counts]");
        h->GetYaxis()->SetTitle("# of entries");
        h_cl_lead_stripq.push_back(h);
        name.str("");
    }

}
//--------------------------------------------------------------//
void VMMAna::drawHistograms()
{
    drawRawHistograms();
    drawClusterMultiplicityHistos();
    drawCleanedChargeHistos();
    drawClusterChargeHistos();
    drawClusterPositionHistos();
    drawClusterChargeHistos_OR();
    drawClusterPositionHistos_OR();
    drawClusterResidualHistos();
    drawCorrHistos();
    //drawChargeRatioHistos();
    //drawLeadStripChargeHistos();
}
//--------------------------------------------------------------//
void VMMAna::drawRawHistograms()
{
    stringstream save_name;
    save_name << "pdo_raw_" << runNumberStr() << ".eps";
    // raw pdo over all channels and vmms
    c_pdo_raw = new TCanvas("c_pdo_raw", "", 1400, 800);
    c_pdo_raw->Divide(4,2);
    for(int i = 0; i < 8; i++) {
        c_pdo_raw->cd(i+1);
        h_pdo_raw.at(i)->Draw();
    }
    c_pdo_raw->SaveAs(save_name.str().c_str()); save_name.str("");

    // raw vmm channel hits
    save_name << "vmm_channel_hits_raw_" << runNumberStr() << ".eps";
    c_vmm_channel_hits_raw = new TCanvas("c_vmm_channel_hits_raw", "", 1000, 800);
    c_vmm_channel_hits_raw->Divide(4,2);
    for(int i = 0; i < 8; i++) {
        c_vmm_channel_hits_raw->cd(i+1);
        h_vmm_hits_raw.at(i)->Draw();
    }
    c_vmm_channel_hits_raw->SaveAs(save_name.str().c_str()); save_name.str("");

    // raw T-chamber hits
    save_name << "strip_hits_raw_" << runNumberStr() << ".eps";
    c_strip_hits_raw = new TCanvas("c_strip_hits_raw", "", 1000, 800);
    c_strip_hits_raw->Divide(1,2);
    for(int i = 0; i < 2; i++) {
        c_strip_hits_raw->cd(i+1);
        h_strip_hits_raw.at(i)->Draw();
    }
    c_strip_hits_raw->SaveAs(save_name.str().c_str());

}
//--------------------------------------------------------------//
void VMMAna::drawClusterMultiplicityHistos()
{
    stringstream save_name;
    c_cluster_mult = new TCanvas("c_cluster_mult", "", 800, 400);
    c_cluster_mult->Divide(3,2);

    for(int i = 0; i < 3; i++) {
        if(i==0) {
            c_cluster_mult->cd(i+1);
            h_num_cluster.at(0)->Draw();
            c_cluster_mult->cd(i+4);
            h_num_cluster.at(1)->Draw();
        }
        else if(i==1) {
            c_cluster_mult->cd(i+1);
            h_cluster_size.at(0)->Draw();
            c_cluster_mult->cd(i+4);
            h_cluster_size.at(1)->Draw();
        }
        else if(i==2) {
            c_cluster_mult->cd(i+1);
            h2_num_cluster_vs_cluster_size.at(0)->Draw("colz");
            c_cluster_mult->cd(i+4);
            h2_num_cluster_vs_cluster_size.at(1)->Draw("colz");
        }
    }

   // c_cluster_mult->cd(1);
   // h_num_cluster->Draw();
   // c_cluster_mult->cd(2);
   // h_cluster_size->Draw();
   // c_cluster_mult->cd(3);
   // h2_num_cluster_vs_cluster_size->Draw("colz");

    save_name << "cluster_mult_" << runNumberStr() << ".eps";
    c_cluster_mult->SaveAs(save_name.str().c_str());
}

//--------------------------------------------------------------//
void VMMAna::drawCleanedChargeHistos()
{
    stringstream save_name;
    c_pdo = new TCanvas("c_pdo", "", 1200, 600);
    c_pdo->Divide(2,2);

    c_pdo->cd(1);
    h_pdo.at(0)->Draw();
    c_pdo->cd(3);
    h_pdo.at(1)->Draw();

    
    c_pdo->cd(2);
    h_pdo_calib.at(0)->Draw();
    c_pdo->cd(4);
    h_pdo_calib.at(1)->Draw();

   // c_pdo->cd(3);
   // h_pdo_adc_fix.at(0)->Draw();
   // c_pdo->cd(6);
   // h_pdo_adc_fix.at(1)->Draw();

    save_name << "strip_charge_" << runNumberStr() << ".eps";
    c_pdo->SaveAs(save_name.str().c_str());

}
//--------------------------------------------------------------//
void VMMAna::drawClusterChargeHistos()
{
    stringstream save_name;

    c_cluster_charge = new TCanvas("c_cluster_charge", "", 1200, 600);
    c_cluster_charge->Divide(2,2);

    ////////////// T6
    // cluster charge
    c_cluster_charge->cd(1);
    h_cl_charge.at(0)->Draw();
    c_cluster_charge->cd(3);
    h_cl_charge.at(1)->Draw();

    // w/ calib
    c_cluster_charge->cd(2);
    h_cl_charge_calib.at(0)->Draw();
    c_cluster_charge->cd(4);
    h_cl_charge_calib.at(1)->Draw();

  //  // with adc smearing
  //  c_cluster_charge->cd(3);
  //  h_cl_charge_adc_fix.at(0)->Draw();
  //  c_cluster_charge->cd(6);
  //  h_cl_charge_adc_fix.at(1)->Draw();

    save_name << "cluster_charge_" << runNumberStr() << ".eps";
    c_cluster_charge->SaveAs(save_name.str().c_str());

}
//--------------------------------------------------------------//
void VMMAna::drawClusterChargeHistos_OR()
{
    stringstream save_name;
    c_cluster_charge_OR = new TCanvas("c_cluster_charge_OR", "", 1200, 600);
    c_cluster_charge_OR->Divide(2,2);

    /////////////// T6
    // charge
    c_cluster_charge_OR->cd(1);
    h_cl_charge_OR.at(0)->Draw();
    c_cluster_charge_OR->cd(3);
    h_cl_charge_OR.at(1)->Draw();

    // w/ claib
    c_cluster_charge_OR->cd(2);
    h_cl_charge_calib_OR.at(0)->Draw();
    c_cluster_charge_OR->cd(4);
    h_cl_charge_calib_OR.at(1)->Draw();

  //  // with adc fix/smearing
  //  c_cluster_charge_OR->cd(3);
  //  h_cl_charge_adc_fix_OR.at(0)->Draw();
  //  c_cluster_charge_OR->cd(6);
  //  h_cl_charge_adc_fix_OR.at(1)->Draw();

    save_name << "cluster_charge_POSTOR_" << runNumberStr() << ".eps";
    c_cluster_charge_OR->SaveAs(save_name.str().c_str());

}
//--------------------------------------------------------------//
void VMMAna::drawClusterPositionHistos()
{
    stringstream save_name;
    c_cluster_position = new TCanvas("c_cluster_position", "", 1200, 600);
    c_cluster_position->Divide(1,2);

    // T6
    c_cluster_position->cd(1);
    h_cl_position.at(0)->Draw();

    // T7
    c_cluster_position->cd(2);
    h_cl_position.at(1)->Draw();

    save_name << "cluster_position_" << runNumberStr() << ".eps";
    c_cluster_position->SaveAs(save_name.str().c_str());
}

//--------------------------------------------------------------//
void VMMAna::drawClusterPositionHistos_OR()
{
    stringstream save_name;
    c_cluster_position_OR = new TCanvas("c_cluster_position_OR", "", 1200, 600);
    c_cluster_position_OR->Divide(1,2);

    // T6
    c_cluster_position_OR->cd(1);
    h_cl_position_OR.at(0)->Draw();

    // T7
    c_cluster_position_OR->cd(2);
    h_cl_position_OR.at(1)->Draw();

    save_name << "cluster_position_POSTOR_" << runNumberStr() << ".eps";
    c_cluster_position_OR->SaveAs(save_name.str().c_str());
}
//--------------------------------------------------------------//
void VMMAna::drawClusterResidualHistos()
{
    stringstream save_name;
    c_cluster_residual = new TCanvas("c_cluster_residual", "", 800,600);
    h_cl_residual->Draw();

    h_cl_residual->SetName("residual");
    h_cl_residual->SaveAs("residual_histo.root");

    save_name << "cluster_position_residual_" << runNumberStr() << ".eps";
    c_cluster_residual->SaveAs(save_name.str().c_str());
    save_name.str("");
    save_name << "cluster_position_residual_" << runNumberStr() << ".root";
    c_cluster_residual->SaveAs(save_name.str().c_str());


    save_name.str("");
    c_t6cl_vs_t7cl = new TCanvas("c_t6cl_vs_t7cl", "", 800, 600);
    c_t6cl_vs_t7cl->cd();
    h2_cl_position_T6vsT7->Draw("colz");
    save_name << "cl_position_T6vsT7_" << runNumberStr() << ".root";
    c_t6cl_vs_t7cl->SaveAs(save_name.str().c_str());
}
//--------------------------------------------------------------//
void VMMAna::drawCorrHistos()
{
    stringstream save_name;
    c_t6vst7_charge = new TCanvas("c_t6vst7_charge", "", 800, 600);
    h2_cl_charge_T6vsT7->Draw("colz"); 

    save_name << "cl_charge_T6vsT7_" << runNumberStr() << ".eps";
    c_t6vst7_charge->SaveAs(save_name.str().c_str());

    // chamber strip correlation
    save_name.str("");
    c_t6vst7_strip = new TCanvas("c_t6vst7_strip", "", 800, 600);
    h2_cl_strip_T6vsT7->Draw("colz");

    save_name << "chamber_strip_T6vsT7_" << runNumberStr() << ".root";
    c_t6vst7_strip->SaveAs(save_name.str().c_str());
}
//--------------------------------------------------------------//
void VMMAna::drawChargeRatioHistos()
{
    stringstream save_name;
    c_cl_strip_qratio = new TCanvas("c_cl_strip_qratio", "", 1200, 800);

    c_cl_strip_qratio->Divide(4,2);

    for(int i = 0; i < (int)h_cl_strip_qratio.size(); i++) {
        c_cl_strip_qratio->cd(i+1);
        h_cl_strip_qratio.at(i)->Draw();
    }

    save_name << "cl_strip_qratio_" << runNumberStr() << ".eps";
    c_cl_strip_qratio->SaveAs(save_name.str().c_str());
}
//--------------------------------------------------------------//
void VMMAna::drawLeadStripChargeHistos()
{
    stringstream save_name;
    c_cl_lead_stripq = new TCanvas("c_cl_lead_stripq", "", 800, 600);
    c_cl_lead_stripq->Divide(2,1);

    for(int i = 0; i < (int)h_cl_lead_stripq.size(); i++) {
        c_cl_lead_stripq->cd(i+1);
        h_cl_lead_stripq.at(i)->Draw();
    } // i

    save_name << "cl_lead_stripQ_" << runNumberStr() << ".eps";
    c_cl_lead_stripq->SaveAs(save_name.str().c_str());

}
//--------------------------------------------------------------//
// load this event's hits
int VMMAna::loadChamberHits()
{
    hits.clear();
    int n_chip = (int)m_chip->size();
    for(int ichip = 0; ichip < n_chip; ichip++) {
        int n_hit_for_chip = (int)m_tdo->at(ichip).size();
        for(int ihit = 0; ihit < n_hit_for_chip; ihit++) {
            vmm::Hit h;
            int chip_ = m_chip->at(ichip);
            int bcid_ = m_bcid->at(ichip).at(ihit);
            int gray_ = m_grayDecoded->at(ichip).at(ihit);
            int channel_ = m_channel->at(ichip).at(ihit);
            int strip_ = TChamberMapping(chip_, channel_);

            // for T7 chamber, rotate the strips about center to align with T6
            if(chip_>=4 && chip_<8) strip_ = 257-strip_;

            int threshold_ = m_threshold->at(ichip).at(ihit);
            int pdo_ = m_pdo->at(ichip).at(ihit);
            int tdo_ = m_tdo->at(ichip).at(ihit);

            if(chip_>=0 && chip_<4) {
                //// fix offset in T6 by +2
                //if(n_T6offset_counter<10) {
                //    cout << "VMMAna::loadChamberHits    Offseting T6 strip number by -2" << endl;
                //    n_T6offset_counter++;
                //}
                //strip_ = strip_ - 2;
                h.setChamberNo(0);
            }
            else if(chip_>=4 && chip_<8) {
                h.setChamberNo(1);
            }

            h.fill(chip_, bcid_, gray_, channel_, strip_, threshold_, pdo_, tdo_);

            hits.push_back(h);
        }
    }
    return (int)hits.size();
}
//--------------------------------------------------------------//
void VMMAna::removeDuplicateHits()
{
    bool removeFirstDuplicate = false;
    int n_duplicates_per_channel = 0;

    bool move_i = true;
    bool move_j = true;
    if(hits.size()>1) {
        for(int i = 0; i < (int)hits.size(); i++ ) {
            if(!move_i) i--;
            move_i = true;

            vmm::Hit hiti = hits.at(i);

            n_duplicates_per_channel = 0;
            for(int j = i+1; j < (int)hits.size(); j++) {
                if(!move_j) j--;
                move_j = true;

                vmm::Hit hitj = hits.at(j);

                if(hiti==hitj) {
                    if(m_dbg) {
                        cout << "VMMAna::removeDuplicateHits    "
                             << "Removing hit[" << j << "] overlapping with hit[" << i << "]"
                             << "  [" << j << "]: (chip, pdo, tdo, bcid, gray) = ("
                             << hitj.chip() << ", " << hitj.pdo() << ", " << hitj.tdo()
                             << ", " << hitj.bcid() << ", " << hitj.gray() << ") "
                             << " [" << i << "]: (" << hiti.chip() << ", " << hiti.pdo()
                             << ", " << hiti.tdo() << ", " << hiti.bcid() << ", " << hiti.gray()
                             << ")" << endl;
                    }

                    hits.erase(hits.begin() + j);
                    move_j = false;
                    //j--;
                    n_duplicates_per_channel++;
                } // duplicate
                else if(hiti.isCrossFire(hitj)) {
                    if(hiti.channel()>hitj.channel()) {
                        if(m_dbg) {
                            cout << "VMMAna::removeDuplicateHits    "
                                 << "Removing CROSSFIRE hit j [" << j << "] overlapping with hit i [" << i << "]"
                                 << "  [" << j << "]: (chip, pdo, tdo, bcid, gray) = ("
                                 << hitj.chip() << ", " << hitj.pdo() << ", " << hitj.tdo()
                                 << ", " << hitj.bcid() << ", " << hitj.gray() << ") "
                                 << " [" << i << "]: (" << hiti.chip() << ", " << hiti.pdo()
                                 << ", " << hiti.tdo() << ", " << hiti.bcid() << ", " << hiti.gray()
                                 << ")" << endl;
                        }

                        hits.erase(hits.begin() + j);
                        move_j = false;
                        //j--;
                    } // i > j channel
                    else {
                        if(m_dbg) {
                            cout << "VMMAna::removeDuplicateHits    "
                                 << "Removing CROSSFIRE hit i [" << i << "] overlapping with hit j [" << j << "]"
                                 << "  [" << i << "]: (chip, pdo, tdo, bcid, gray) = ("
                                 << hiti.chip() << ", " << hiti.pdo() << ", " << hiti.tdo()
                                 << ", " << hiti.bcid() << ", " << hiti.gray() << ") "
                                 << " [" << j << "]: (" << hitj.chip() << ", " << hitj.pdo()
                                 << ", " << hitj.tdo() << ", " << hitj.bcid() << ", " << hitj.gray()
                                 << ")" << endl;
                        }
                        hits.erase(hits.begin() + i);
                        //i--;
                        move_i = false;
                    }
                } // cross-fire
            } // j
        } // i
    }
}

//--------------------------------------------------------------//
// build up clusters from all of the hits
int VMMAna::makeClusters()
{

    std::vector<vmm::Hit> tmpHits;
    for(int i = 0; i < (int)hits.size(); i++) {
        tmpHits.push_back(hits.at(i));
    }

    bool move_j = true;
    clusters.clear();
    for(int i = 0; i < (int)tmpHits.size(); i++) {
        vmm::Cluster cluster;
        vmm::Hit hiti = tmpHits.at(i);

        int chamber_no = -1;
        if(hiti.chip()>=0 && hiti.chip()<4) {
            chamber_no = 0; // T6
        }
        else if(hiti.chip()>=4 && hiti.chip()<8) {
            chamber_no = 1; // T7
        }

        if(chamber_no < 0) {
            cout << "VMMAna::makeClusters    Cluster at negative chamber number. Exitting." << endl;
            exit(1);
        }

        cluster.setChamberNo(chamber_no);
        cluster.addHit(hiti);

        for(int j = i+1; j < (int)tmpHits.size(); j++) {
            if(!move_j) j--;
            move_j = true;

            vmm::Hit hitj = tmpHits.at(j);
            int j_chamber = -1;
            if( hitj.chip()>=0 && hitj.chip()<4) j_chamber = 0;
            else if(hitj.chip()>=4 && hitj.chip()<8) j_chamber = 1;
            if(j_chamber<0) {
                cout << "VMMAna::makeClusters Cluster j at negative chamber number. Exitting." << endl;
                exit(1);
            }
            if( !(chamber_no==j_chamber) ) continue;


            if( (fabs(hiti.gray()-hitj.gray())<=12) &&
                (fabs(hiti.strip()-hitj.strip())<=5) &&
                (fabs(hiti.strip()-hitj.strip())>0) ) {

                cluster.addHit(hitj);
                tmpHits.erase(tmpHits.begin() + j);
                move_j = false;
                //j--;
            } // match
            
        } // j
        clusters.push_back(cluster);
    } // i



 //   for(int i = (int)tmpHits.size()-1; i >= 0; i--) {
 //       vmm::Cluster cluster;
 //       vmm::Hit hiti = tmpHits.at(i);
 //       cluster.addHit(hiti);


 //       for(int j = i-1; j>=0; j--) {
 //           vmm::Hit hitj = tmpHits.at(j);

 //           if( (fabs(hiti.gray()-hitj.gray())<=12) &&
 //               (fabs(hiti.strip()-hitj.strip())<=5) &&
 //               (fabs(hiti.strip()-hitj.strip())>0) ) {

 //               cluster.addHit(hitj);

 //               tmpHits.erase(tmpHits.begin() + j);
 //           }
 //       } // j

 //       // ith cluster
 //       clusters.push_back(cluster);
 //   } // i

    return clusters.size();
}
//--------------------------------------------------------------//
void VMMAna::getClusterChargeAndPosition()
{
    if(clusters.size()>0) {
        for(auto& cluster : clusters) {

            int cl_pdo = cluster.pdo();
            double cl_pdo_calibrated = 0.0;
            double cl_position_centroid = 0.0;
            double cl_position_centroid_calibrated = 0.0;

            std::vector<vmm::Hit>& cluster_strips = cluster.hits();

            for(auto strip : cluster_strips) {
                double strip_pdo_calibrated;
                if(calib()) {
                    cout << "VMMAna::fillClusterHistos    Attempting to apply calibration. This is"
                        << " not setup. Exitting." << endl;
                    exit(1);
                }
                else {
                    strip_pdo_calibrated = (GAIN*strip.pdo()) / m_channel_gain_calib[cluster.chamber()][strip.channel()];
                }

                // CALIBRATED PDO
                cl_pdo_calibrated += strip_pdo_calibrated; 

                // CENTROID
                cl_position_centroid += strip.pdo() * (strip.strip()*PITCH + PITCH/2.);
                cl_position_centroid_calibrated += strip_pdo_calibrated * (strip.strip()*PITCH + PITCH/2.);

            } // strip

            cluster.setPdoCalibrated(cl_pdo_calibrated);
            cluster.setPositionCentroid(cl_position_centroid / cl_pdo);
            cluster.setPositionCentroidCalibrated(cl_position_centroid_calibrated / cl_pdo_calibrated);

        } // cluster

    } // clusters > 0


}
//--------------------------------------------------------------//
void VMMAna::removeDuplicateClusters()
{

 //   int initial_size = (int)clusters.size();
 //   cout << "-----------------------------------" << endl;
    std::sort(clusters.begin(), clusters.end(), byClusPdo);
 //   cout << initial_size << " clusters before removal: " << endl;
 //   for(int i = 0; i < (int)clusters.size(); i++) clusters.at(i).print();
    //for(int i = 0; i < (int)clusters.size(); i++) {
    //    cout << clusters.at(i).pdo() << " ";
    //}
    //cout << endl;
    clusters.erase(std::unique(clusters.begin(), clusters.end(), clusPdoEqual), clusters.end());

  //  int final_size = (int)clusters.size();   
 
    //cout << "clusters AFTER removal: ";
    //for(int i = 0; i < (int)clusters.size(); i++) {
    //    cout << clusters.at(i).pdo() << " ";
    //}
    //cout << endl;
  //  if(initial_size>final_size) {
  //      cout << "CLUSTER REMOVED -- final size: " << final_size << endl;
  //      for(int i = 0; i < (int)clusters.size(); i++) clusters.at(i).print();
  //  }

 //   static std::set<vmm::Cluster> clustersToRemove
 //   clustersToRemove.clear();

 //   for(int i = 0; i < (int)clusters.size(); ) {
 //       vmm::Cluster clusteri = clusters.at(i);
 //       for(int j=i+1; j < (int)clusters.size(); j++) {

 //           vmm::Cluster clusterj = clusters.at(j);

 //           if(clusteri.pdo() == clusterj.pdo()) clustersToRemove.insert(clusterj);

 //       } // j

 //   } // i
//std::unique (myvect.begin(), myvect.end(),
//             [](const myclass& a, const myclass& b) {
//               return /* some expression to test equality */;
//             }); 
}
//--------------------------------------------------------------//
void VMMAna::removeDuplicateStripsInClusters()
{
    if(clusters.size()>0) {
        for(auto& cluster : clusters) {
            cluster.sortClusterHitsByChannel();
            cluster.removeDuplicateStrips();
        } // cluster

    }

}
//--------------------------------------------------------------//
// TSelector::Terminate
void VMMAna::Terminate()
{
    cout << " --------------------- VMMAna::Terminate ---------------------- " << endl;

    drawHistograms();

    cout << "Total # of events passing all cuts: " << n_total << endl;
}

//--------------------------------------------------------------//
void VMMAna::initializeNtupleVariablesAndBranches()
{
    if(!m_chain) {
        cout << " ------ VMMAna::initializeNtupleVariablesAndBranches ------ " << endl;
        cout << "   Input chain is null. You must set this using the method  " << endl;
        cout << "   \"VMMAna::setChain(TChain* chain)\" method at the start  " << endl;
        cout << "   of your job. >>> Exiting.                                " << endl;
        cout << " ---------------------------------------------------------- " << endl;
        exit(1);
    }

    // clear
    m_eventFAFA         = 0;
    m_triggerTimeStamp  = 0;
    m_triggerCounter    = 0;
    m_chip              = 0;
    m_eventSize         = 0;
    m_tdo               = 0;
    m_pdo               = 0;
    m_flag              = 0;
    m_threshold         = 0;
    m_bcid              = 0;
    m_grayDecoded       = 0;
    m_channel           = 0;

    m_chain->SetBranchAddress("eventFAFA", &m_eventFAFA, &b_eventFAFA);
    m_chain->SetBranchAddress("triggerTimeStamp", &m_triggerTimeStamp, &b_triggerTimeStamp);
    m_chain->SetBranchAddress("triggerCounter", &m_triggerCounter, &b_triggerCounter);
    m_chain->SetBranchAddress("chip", &m_chip, &b_chip);
    m_chain->SetBranchAddress("eventSize", &m_eventSize, &b_eventSize);
    m_chain->SetBranchAddress("tdo", &m_tdo, &b_tdo);
    m_chain->SetBranchAddress("pdo", &m_pdo, &b_pdo);
    m_chain->SetBranchAddress("flag", &m_flag, &b_flag);
    m_chain->SetBranchAddress("threshold", &m_threshold, &b_threshold);
    m_chain->SetBranchAddress("bcid", &m_bcid, &b_bcid);
    m_chain->SetBranchAddress("grayDecoded", &m_grayDecoded, &b_grayDecoded);
    m_chain->SetBranchAddress("channel", &m_channel, &b_channel);

}
//--------------------------------------------------------------//
// mapping
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
