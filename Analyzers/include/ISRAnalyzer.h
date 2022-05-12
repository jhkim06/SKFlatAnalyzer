#ifndef ISRAnalyzer_h
#define ISRAnalyzer_h

#include "TUnfoldBinning.h"
#include "ISRUnfold.h"
#include "SMPAnalyzerCore.h"

enum class DressedMode{
    /*
     AllPHOTON : add all stable photons
     MotherMatch : add stable photons matched to DY history
     DRMatch : add stable photons matched with DR cut
     
     Study the DR distribution when MotherMatch used and use it as delta R matching condition
     */
    AllPHOTON = 0, MotherMatch, DRMatch, MotherDRMatch
};

enum GenMode{
    PostFSR = 0, DressedDRp1, DressedDRp4, PreFSR
};

enum GenPID{
    TOP = 6, ELECTRON = 11, MUON = 13, TAU = 15, PHOTON = 22, PROTON = 2212
};

const int nmass_window = 6;
const double mass_window[] = {40, 64, 81, 101, 200, 320, 1000};
const double mass_window_an026[] = {40, 76, 106, 170, 350, 1000};

const int n_pt_bin_fine=36;
const double pt_bin_fine[]={0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 16., 18., 20.5, 23, 25.5, 28., 31., 34., 37., 40., 43.75, 47.5, 51.25, 55., 60., 65., 70., 75., 81.25, 87.5, 93.75, 100.};
                                                                                                       
//const double pt_bin[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 25, 28, 32, 37, 43, 52, 65, 80, 100, 130, 160, 190, 220, 250, 300, 350, 400, 450, 500, 1000}; // 27 bins
const double pt_bin[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 25, 28, 32, 37, 43, 52, 65, 80, 100}; // 27 bins
//const double ptbin[]={0, 2, 4, 6, 8, 10, 12, 14, 18, 22, 28, 37, 52, 65, 100}; // 14 bins
//const double ptbin[]={0, 2, 4, 6, 8, 10, 12, 14, 18, 22, 28, 52, 100}; // 12 bins

const int n_pt_bin_coarse=18;
const double pt_bin_coarse[]={0., 2., 4., 6., 8., 10., 12., 14., 18., 23, 28., 34., 40., 47.5, 55., 65., 75., 87.5, 100.};

const int n_mass_bin_fine_muon = 72;
const double mass_bin_fine_muon[n_mass_bin_fine_muon+1] = {40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000};

const int n_mass_bin_fine_electron = 68;
const double mass_bin_fine_electron[n_mass_bin_fine_electron+1] = {50, 52.5, 55, 57.5,60, 62,64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000};

const int n_mass_bin_coarse_muon= 36;
const double mass_bin_coarse_muon[n_mass_bin_coarse_muon+1] = {40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000};

const int n_mass_bin_coarse_electron = 34;
const double mass_bin_coarse_electron[n_mass_bin_coarse_electron+1] = {50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000};

class ISRAnalyzer : public ISRUnfold  {

public:

    void initializeAnalyzer();
    void executeEventFromParameter(AnalyzerParameter param);
    void executeEventWithChannelName(TString channelname);
    void executeEvent();

    void fill_ISR_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> map_weight);

    //
    int get_DY_gen_particles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& letpon0, Gen& lepton1, int mode);
    int get_DY_gen_particles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& letpon0, Gen& lepton1, int mode, vector<const Gen*>& added_photons);
    int get_DY_bare_lepton_pair(const vector<Gen>& gens, const vector<const Gen*>& leptons, Gen& lepton0, Gen& lepton1, bool verbose=false);
    int get_DY_dressed_lepton_pair(const vector<Gen>& gens, const vector<const Gen*>& leptons, vector<const Gen*>& photons, Gen& lepton0, Gen& lepton1,
                                    const DressedMode mode, vector<const Gen*>& added_photons, const double dR = 0.1);
    void save_gen_history(const vector<Gen>& gens, const Gen& lepton, vector<int>& partindex_vector, const int index_limit = -1);
    void print_gen_particles(const vector<Gen>& gens);
    
    ISRAnalyzer();
    ~ISRAnalyzer();
    
    class ISRParameter : public Parameter{
        
    public:
        double dilep_pt_cut = 100;
        
        inline ISRParameter Clone(vector<Lepton*> leps_,int weightbit_=-1){
            ISRParameter out=*this;
            out.leps=leps_;
            if(weightbit_>=0) out.weightbit=weightbit_;
            return out;
        }

        ISRParameter(){
            Parameter();
        }
        ISRParameter(TString elID, vector<TString> Trig, double l0ptcut=-1, double l1ptcut=-1, double diptcut = 100, vector<Lepton*> leps_={}){
            Parameter(elID, Trig, l0ptcut, l1ptcut, leps_);
            dilep_pt_cut = diptcut;
        }
        ISRParameter(TString muID, TString muISO, vector<TString> Trig, double l0ptcut=-1, double l1ptcut=-1, double diptcut = 100, vector<Lepton*> leps_={}){
            Parameter(muID, muISO, Trig, l0ptcut, l1ptcut, leps);
            dilep_pt_cut = diptcut;
        }
    };
    
private:
    bool IsNominalRun=true;
    bool IsSkimmed=false;
    
    TUnfoldParameter* tunfold_parameter;
    //ISRUnfold tunfold_hists;
    //ISRUnfold* tunfold_hists; //
    
    // tunfold_pt_2d, tunfold_mass_2d, tunfold_pt_1d, tunfold_mass_1d

};




#endif

