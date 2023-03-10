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
    AllPHOTON=0, MotherMatch, DRMatch, MotherDRMatch
};

enum GenMode{
    PostFSR=0, DressedDRp1, DressedDRp4, PreFSR
};

enum GenPID{
    TOP=6, ELECTRON=11, MUON=13, TAU=15, PHOTON=22, PROTON=2212
};

const int nmass_window = 6;
const double mass_window[] = {40, 64, 81, 101, 200, 320, 1000};
const double mass_window_an026[] = {40, 76, 106, 170, 350, 1000};

                                                                                                       
const double pt_bin_extended[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 25, 28, 32, 37, 43, 52, 65, 80, 100, 130, 160, 190, 220, 250, 300, 350, 400, 450, 500, 1000}; // 38 bins
const double pt_bin[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 25, 28, 32, 37, 43, 52, 65, 80, 100}; // 27 bins
//const double ptbin[]={0, 2, 4, 6, 8, 10, 12, 14, 18, 22, 28, 37, 52, 65, 100}; // 14 bins
//const double ptbin[]={0, 2, 4, 6, 8, 10, 12, 14, 18, 22, 28, 52, 100}; // 12 bins

const int n_pt_bin_fine=36;
const double pt_bin_fine[]={0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 16., 18., 20.5, 23, 25.5, 28., 31., 34., 37., 40., 43.75, 47.5, 51.25, 55., 60., 65., 70., 75., 81.25, 87.5, 93.75, 100.};

const int n_pt_bin_coarse=18;
const double pt_bin_coarse[]={0., 2., 4., 6., 8., 10., 12., 14., 18., 23, 28., 34., 40., 47.5, 55., 65., 75., 87.5, 100.};

const int n_pt_bin_fine_v2=30;
const double pt_bin_fine_v2[]={0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.5, 23, 25.5, 28., 31., 34., 37., 40., 43.75, 47.5, 51.25, 55., 60., 65., 70., 75., 81.25, 87.5, 93.75, 100.};

const int n_pt_bin_coarse_v2=15;
const double pt_bin_coarse_v2[]={0., 4., 8., 12., 16., 18., 23, 28., 34., 40., 47.5, 55., 65., 75., 87.5, 100.};

const int n_pt_bin_fine_v3=30;
const double pt_bin_fine_v3[]={0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.5, 23, 25.5, 28., 31., 34., 37., 40., 43.75, 47.5, 51.25, 55., 60., 65., 70., 75., 81.25, 87.5, 93.75, 100.};

const int n_pt_bin_coarse_v3=10;
const double pt_bin_coarse_v3[]={0., 4., 8., 12., 16., 18., 28., 40., 55., 75., 100.};

const int n_pt_bin_fine_v4=40;
const double pt_bin_fine_v4[]={0., 2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5, 30., 32.5, 35., 37.5, 40., 42.5, 45., 47.5, 50., 52.5, 55., 57.5, 60., 62.5, 65., 67.5, 70., 72.5, 75., 77.5, 80., 82.5, 85., 87.5, 90., 92.5, 95., 97.5, 100};

const int n_pt_bin_coarse_v4=20;
const double pt_bin_coarse_v4[]={0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100};

const int n_mass_bin_fine_mu = 72;
const double mass_bin_fine_mu[n_mass_bin_fine_mu+1] = {40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000};

const int n_mass_bin_fine_el = 68;
const double mass_bin_fine_el[n_mass_bin_fine_el+1] = {50, 52.5, 55, 57.5,60, 62,64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000};

const int n_mass_bin_coarse_mu= 36;
const double mass_bin_coarse_mu[n_mass_bin_coarse_mu+1] = {40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000};

const int n_mass_bin_coarse_el = 34;
const double mass_bin_coarse_el[n_mass_bin_coarse_el+1] = {50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000};

class ISRAnalyzer : public ISRUnfold  {
public:

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
        // electron
        ISRParameter(TString elID, vector<TString> Trig, double l0ptcut=-1, double l1ptcut=-1, double diptcut = 100, vector<Lepton*> leps_={}){
            Parameter(elID, Trig, l0ptcut, l1ptcut, leps_);
            dilep_pt_cut = diptcut;
        }
        // muon
        ISRParameter(TString muID, TString muISO, vector<TString> Trig, double l0ptcut=-1, double l1ptcut=-1, double diptcut = 100, vector<Lepton*> leps_={}){
            Parameter(muID, muISO, Trig, l0ptcut, l1ptcut, leps);
            dilep_pt_cut = diptcut;
        }
    };

    void initializeAnalyzer();
    void executeEventFromParameter(AnalyzerParameter param);
    void executeEventWithChannelName(TString channelname);
    void executeEvent();

    bool pass_lepton_kinematic_selections(TString channelname, Particle* l0, Particle* l1, const ISRParameter& p);

    void fill_mass_dependent_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> map_weight);
    void fill_dipt_resol_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, Particle* truth_l0, Particle* truth_l1);

    void fill_unfold_isr_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> map_weight, TUnfold_Bin bin_type = TUnfold_Bin::smeared_bin);
    void fill_unfold_isr_responseM(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, Particle* truth_l0, Particle* truth_l1,map<TString,double>& reco_weights, map<TString,double>& gen_weights);
    void fill_unfold_isr_fake_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, Particle* truth_l0, Particle* truth_l1, map<TString,double>& map_weights, const ISRParameter& p);
    void fill_mass_dependent_stability_purity(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, Particle* truth_l0, Particle* truth_l1, map<TString,double> map_weight);

    //
    int get_DY_gen_particles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& letpon0, Gen& lepton1, int mode);
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
    
private:
    bool IsNominalRun=true;
    bool IsSkimmed=false;
    int job_number;
    
    TUnfoldParameter* tunfold_parameter;
    TUnfoldParameter* tunfold_parameter_1d_pt;
    TUnfoldParameter* tunfold_parameter_1d_pt_v2;
    TUnfoldParameter* tunfold_parameter_1d_pt_v3;
    TUnfoldParameter* tunfold_parameter_1d_pt_v4;
    TUnfoldParameter* tunfold_parameter_1d_mass_el;
    TUnfoldParameter* tunfold_parameter_1d_mass_mu;

    map<TString, TUnfoldParameter*> tunfold_parameters;
};

#endif
