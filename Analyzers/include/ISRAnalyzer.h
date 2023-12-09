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

class ISRAnalyzer : public ISRUnfold {
    
public:
    virtual void initializeAnalyzer();
    virtual void executeEvent();
    virtual void executeEventGen(Parameter& p);
    virtual void executeEventWithParameter(Parameter& p);
    virtual void executeEventWithParameter(Parameter&& p){Parameter pp=p;
        executeEventWithParameter(pp);}
    virtual Parameter MakeParameter(TString key, TString option="");
    virtual bool PassSelection(Parameter& p);
    virtual void EvalWeights(Parameter& p);
    virtual void ResetRecoWeights(Parameter& p);
    virtual void FillHists(Parameter& p);
    void FillLeptonHists(Parameter& p);   
    TString add_postfix_BE_EE(Parameter& p, int lepton_index, bool add_sign=false);
    bool PassHEMElectron(Electron electron);
    bool RemoveElectron(Electron electron);
    bool RemoveMuon(Muon muon, double phi_min, double phi_max, bool for_endcap, int eta_sign=0);
    
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

