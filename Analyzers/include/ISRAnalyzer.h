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

const int nmass_window = 6;
const double mass_window_mm[] = {40, 64, 81, 101, 200, 320, 1000};
const double mass_window_ee[] = {50, 64, 81, 101, 200, 320, 1000};
const double mass_window_an026[] = {40, 76, 106, 170, 350, 1000};

const int n_pt_bin_fine=18;
const double pt_bin_fine[]={0., 2., 4., 6., 8., 10., 12., 14., 18., 23, 28., 34., 40., 47.5, 55., 65., 75., 87.5, 100.};

const int n_pt_bin_coarse=9;
const double pt_bin_coarse[]={0., 4., 8., 12., 18., 28., 40., 55., 75., 100.};

const int n_mass_bin_fine_mu = 72;
const double mass_bin_fine_mu[n_mass_bin_fine_mu+1] = {40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000};

const int n_mass_bin_fine_el = 68;
const double mass_bin_fine_el[n_mass_bin_fine_el+1] = {50, 52.5, 55, 57.5,60, 62,64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000};

const int n_mass_bin_coarse_mu= 36;
const double mass_bin_coarse_mu[n_mass_bin_coarse_mu+1] = {40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000};

const int n_mass_bin_coarse_el = 34;
const double mass_bin_coarse_el[n_mass_bin_coarse_el+1] = {50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000};

class ISRAnalyzer : public ISRUnfold {

public:
  virtual void initializeAnalyzer();
  virtual void executeEvent();
  virtual void executeEventGen();
  virtual void executeEventWithParameter(Parameter& p);
  virtual void executeEventWithParameter(Parameter&& p){Parameter pp=p;executeEventWithParameter(pp);}
  virtual Parameter MakeParameter(TString key);
  virtual bool PassSelection(Parameter& p);
  virtual void EvalWeights(Parameter& p);
  virtual void ResetRecoWeights(Parameter& p);
  virtual void FillHists(Parameter& p);

    bool pass_lepton_kinematic_selections(const Parameter& p, Particle* l0, Particle* l1);
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

  //void FillGenAFBHists(TString pre,TString suf,const Gen& genl0,const Gen& genl1,const Gen& genphotons,double w);

  TString hardprefix;
  map<TString,TH3D*> map_hist_cost;
  bool IsNominalRun=true;
  bool IsSkimmed=false;

private:
  int job_number;

  TUnfoldParameter* tunfold_parameter_ee;
  TUnfoldParameter* tunfold_parameter_mm;

  map<TString, TUnfoldParameter*> tunfold_parameters;

};



#endif

