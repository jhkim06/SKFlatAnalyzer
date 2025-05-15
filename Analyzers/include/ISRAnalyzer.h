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
    virtual Parameter MakeParameter(TString key, TString option="", TString option_shared="");
    std::vector<Electron> ISRGetElectrons(TString id, double ptmin, double fetamax, bool vetoHEM = false, int set=0, int member=0) const; 
    virtual bool PassSelection(Parameter& p);
    virtual void EvalWeights(Parameter& p);
    virtual void ResetRecoWeights(Parameter& p);
    virtual void FillHists(Parameter& p);
    void FillLeptonHists(Parameter& p);   
    TString add_postfix_BE_EE(Parameter& p, int lepton_index, bool add_sign=false);
    bool PassHEMElectron(Electron electron);
    bool RemoveElectron(Electron electron);
    bool RemoveMuon(Muon muon, double phi_min, double phi_max, bool for_endcap, int eta_sign=0);
    
    void get_DY_gen_particles();
    void get_DY_bare_lepton_pair(bool verbose=false);
    void get_DY_dressed_lepton_pair(Gen& lepton0, Gen& lepton1, const DressedMode mode, const double dR = 0.1);
    void save_gen_history(const vector<Gen>& gens, const Gen& lepton, vector<int>& partindex_vector, const int index_limit = -1);
    void print_gen_particles(const vector<Gen>& gens);
    
    void set_job_number(int job_num){
        job_number=job_num;
    }
    
    ISRAnalyzer();
    ~ISRAnalyzer();
    
    TString hardprefix;
    bool IsNominalRun=true;
    bool IsSkimmed=false;

    Gen gen_p0_isr, gen_p1_isr, gen_l0_bare, gen_l1_bare;
    vector<Gen*> gen_leptons;
    vector<Gen*> gen_photons;
    int DY_index;
    
private:
    int job_number;
};

#endif

