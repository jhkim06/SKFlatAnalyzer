#ifndef ISRAnalyzer_h
#define ISRAnalyzer_h

#include "TKey.h"
#include "SMPAnalyzerCore.h"
#include "TUnfoldBinning.h"
#include "ISRUnfold.h"

enum class DressedMode{
    /*
     AllPHOTON : add all stable photons
     MotherMatch : add stable photons matched to DY history
     DRMatch : add stable photons matched with DR cut
     
     Study the DR distribution when MotherMatch used and use it as delta R matching condition
     */
    AllPHOTON=0, MotherMatch, DRMatch, MotherDRMatch
};

enum GenMode{
    PostFSR=0, DressedDRp1, DressedDRp4, PreFSR
};

enum GenPID{
    TOP=6, ELECTRON=11, MUON=13, TAU=15, PHOTON=22, PROTON=2212
};

const vector<double> mass_window_5 = {55., 64., 81., 101., 200., 1000.};
const vector<double> mass_window_5_v1 = {60., 81., 101., 160., 1000.}; // electron
const vector<double> mass_window_5_v2 = {55., 64., 81., 101., 160., 1000.}; // muon 
const vector<double> mass_window_5_v3 = {55., 64., 81., 101., 1000.}; // muon
const vector<double> mass_window_5_v4 = {64., 81., 101., 160., 320.};

// bins for dilepton mass
const vector<double> mass_bin_fine_mu = {40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000};
const vector<double> mass_bin_fine_el = {50, 52.5, 55, 57.5,60, 62,64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000};

const vector<double> mass_bin_coarse_mu = {40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000};
const vector<double> mass_bin_coarse_el = {50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000};


const vector<double> mass_bin_fine = {55, 57.5, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000};
const vector<double> mass_bin_coarse = {55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000};

class ISRAnalyzer : public ISRUnfold {
    
public:
    virtual void initializeAnalyzer();
    virtual void executeEvent();
    virtual void executeEventGen();
    virtual void executeEventWithParameter(Parameter& p);
    virtual void executeEventWithParameter(Parameter&& p){Parameter pp=p;
        executeEventWithParameter(pp);}
    virtual Parameter MakeParameter(TString key);
    virtual bool PassSelection(Parameter& p);
    virtual void EvalWeights(Parameter& p);
    virtual void ResetRecoWeights(Parameter& p);
    virtual void FillHists(Parameter& p);
    
    bool pass_lepton_kinematic_selections(const Parameter& p, Particle* l0, Particle* l1);
    int get_DY_gen_particles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& letpon0, Gen& lepton1);
    int get_DY_gen_particles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& letpon0, Gen& lepton1, int mode, vector<const Gen*>& added_photons);
    int get_DY_bare_lepton_pair(const vector<Gen>& gens, const vector<const Gen*>& leptons, Gen& lepton0, Gen& lepton1, bool verbose=false);
    int get_DY_dressed_lepton_pair(const vector<Gen>& gens, const vector<const Gen*>& leptons, vector<const Gen*>& photons, Gen& lepton0, Gen& lepton1,
                                   const DressedMode mode, vector<const Gen*>& added_photons, const double dR = 0.1);
    void save_gen_history(const vector<Gen>& gens, const Gen& lepton, vector<int>& partindex_vector, const int index_limit = -1);
    void print_gen_particles(const vector<Gen>& gens);
    
    void set_job_number(int job_num){
        job_number=job_num;
    }
    
    ISRAnalyzer();
    ~ISRAnalyzer();
    
    //void FillGenAFBHists(TString pre,TString suf,const Gen& genl0,const Gen& genl1,const Gen& genphotons,double w);
    
    TString hardprefix;
    map<TString,TH3D*> map_hist_cost;
    bool IsNominalRun=true;
    bool IsSkimmed=false;
    
private:
    int job_number;
};

#endif

