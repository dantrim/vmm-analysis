#ifndef T6T7_ANA_H
#define T6T7_ANA_H


//ROOT
#include "TROOT.h"
#include "TChain.h"
#include "TBranch.h"
#include "TFile.h"
#include "TSelector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

//std/stl
#include <vector>
#include <sstream>

//vmmana
#include "hit.hh"
#include "cluster.hh"

namespace t6t7 {

class VMMAna : public TSelector {

    public :
        VMMAna();
        virtual ~VMMAna(){};

        void setChain(TChain* chain) { m_chain = chain; }
        TChain& chain() { return *m_chain; }

        void setDebug(bool dbg) { m_dbg = dbg; }
        bool dbg() { return m_dbg; }

        void setRunNumber(int number) { m_run_number = number; }
        int runNumber() { return m_run_number; }
        std::string runNumberStr() {
            std::stringstream ss;
            ss << std::setw(4) << std::setfill('0') << runNumber();
            return ss.str();
        }

        void setCalibFile(std::string file) { m_calib_file = file; }
        std::string calibFile() { return m_calib_file; }
        void setDoCalib(bool docalib) { m_calib = docalib; }
        bool calib() { return m_calib; }

        void loadCalibration();
        std::vector<std::vector<double> > m_channel_calib;
        std::vector<std::vector<double> > m_channel_gain_calib;
        std::vector<std::vector<double> > m_channel_pedestal_calib;

        // T-Chamber mapping
        int TChamberMapping(int vmm_id, int vmm_channel);

        int loadChamberHits();
        std::vector<vmm::Hit> hits;

        void removeDuplicateHits();

        int makeClusters();
        std::vector<vmm::Cluster> clusters;
        void removeDuplicateClusters();
        void removeDuplicateStripsInClusters();

        std::string intToBinaryStr(int num);
        int binaryStrToInt(std::string binary);

        ///////////////////////////////////////////////
        // initialize all histograms
        ///////////////////////////////////////////////
        void initializeHistograms();
        void initRawHistos();
        void initClusterMultiplicityHistos();
        void initCleanedChargeHistos();
        void initClusterChargeHistos();
        void initClusterPositionHistos();

        ///////////////////////////////////////////////
        // histogram filling
        ///////////////////////////////////////////////
        void fillRawHistograms();
        void fillClusterMultiplicityHistos();
        void fillClusterHistos();

        ///////////////////////////////////////////////
        // histogram drawing
        ///////////////////////////////////////////////
        void drawHistograms();
        void drawRawHistograms();
        void drawClusterMultiplicityHistos();
        void drawCleanedChargeHistos();
        void drawClusterChargeHistos();
        void drawClusterPositionHistos();


        //////////////////////////////////////////////////////////////////////
        // histograms/canvases
        //////////////////////////////////////////////////////////////////////

        // raw histos
        TCanvas* c_pdo_raw;
        std::vector<TH1F*> h_pdo_raw;

        // raw vmm channel hits
        TCanvas* c_vmm_channel_hits_raw;
        std::vector<TH1F*> h_vmm_hits_raw;

        // raw T-chamber strip hits
        TCanvas* c_strip_hits_raw;
        std::vector<TH1F*> h_strip_hits_raw;

        // cluster multiplicity
        TCanvas* c_cluster_mult;
        TH1F* h_num_cluster;
        TH1F* h_cluster_size;
        TH2F* h2_num_cluster_vs_cluster_size;

        // "cleaned" charge histos
        TCanvas* c_pdo;
        std::vector<TH1F*> h_pdo;
        std::vector<TH1F*> h_pdo_adc_fix;
        std::vector<TH1F*> h_pdo_calib;

        // cluster charge histos
        TCanvas* c_cluster_charge;
        std::vector<TH1F*> h_cl_charge;
        std::vector<TH1F*> h_cl_charge_calib;
        std::vector<TH1F*> h_cl_charge_adc_fix;

        // cluster position histos
        TCanvas* c_cluster_position;
        std::vector<TH1F*> h_cl_position;
        








        //////////////////////////////////////////////////////////////////////
        // ---------------------------------------------------------------- //
        //////////////////////////////////////////////////////////////////////
        // ---------------------------------------------------------------- //
        //////////////////////////////////////////////////////////////////////


        ///////////////////////////////////////////////
        // TSelector Methods
        ///////////////////////////////////////////////
        virtual Int_t Version() const { return 2; }
        virtual void Init(TTree* tree);
        virtual Bool_t Notify();
        virtual void Begin(TTree* tree);
        virtual void SlaveBegin(TTree* tree);
        virtual void Terminate();
        virtual Bool_t Process(Long64_t entry);

        ///////////////////////////////////////////////
        // VMM Tree leaves and branches
        ///////////////////////////////////////////////

        // leaves
        int                   m_eventFAFA;
        std::vector<int>           *m_triggerTimeStamp;
        std::vector<int>           *m_triggerCounter;
        std::vector<int>           *m_chip;
        std::vector<int>           *m_eventSize;
        std::vector<std::vector<int> >  *m_tdo;
        std::vector<std::vector<int> >  *m_pdo;
        std::vector<std::vector<int> >  *m_flag;
        std::vector<std::vector<int> >  *m_threshold;
        std::vector<std::vector<int> >  *m_bcid;
        std::vector<std::vector<int> >  *m_grayDecoded;
        std::vector<std::vector<int> >  *m_channel;

        // branches
        TBranch*                b_eventFAFA; //!
        TBranch*                b_triggerTimeStamp; //!
        TBranch*                b_triggerCounter; //!
        TBranch*                b_chip; //!
        TBranch*                b_eventSize; //!
        TBranch*                b_tdo; //!
        TBranch*                b_pdo; //!
        TBranch*                b_flag; //!
        TBranch*                b_threshold; //!
        TBranch*                b_bcid; //!
        TBranch*                b_grayDecoded; //!
        TBranch*                b_channel; //!
        

    private :

        TChain* m_chain; // chain that we are processing
        bool m_dbg;
        std::string m_calib_file;
        bool m_calib;
        int m_run_number; // run number associated with file we are processing 

        ///////////////////////////////////////////////
        // initialization/setup of reading
        ///////////////////////////////////////////////
        void initializeNtupleVariablesAndBranches();
        


}; // class 

} // namespace t6t7

#endif
