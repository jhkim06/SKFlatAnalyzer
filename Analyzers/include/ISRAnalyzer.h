#ifndef ISRAnalyzer_h
#define ISRAnalyzer_h

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

class ISRAnalyzer : public SMPAnalyzerCore  {

public:

    void initializeAnalyzer();
    void executeEventFromParameter(AnalyzerParameter param);
    void executeEvent();

    void FillHists(const TString channelname, const TString pre, const TString suf,
                   const Particle* l0, const Particle* l1, const double weight) const;

    //
    int get_DY_gen_particles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& letpon0, Gen& lepton1, int mode);
    int get_DY_gen_particles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& letpon0, Gen& lepton1, int mode, vector<const Gen*>& added_photons);
    int get_DY_bare_lepton_pair(const vector<Gen>& gens, const vector<const Gen*>& leptons, Gen& lepton0, Gen& lepton1);
    int get_DY_dressed_lepton_pair(const vector<Gen>& gens, const vector<const Gen*>& leptons, vector<const Gen*>& photons, Gen& lepton0, Gen& lepton1,
                                    const DressedMode mode, vector<const Gen*>& added_photons, const double dR = 0.1);
    void save_gen_history(const vector<Gen>& gens, const Gen& lepton, vector<int>& partindex_vector, const int index_limit = -1);
    void print_gen_particles(const vector<Gen>& gens);
    
    ISRAnalyzer();
    ~ISRAnalyzer();

};

enum GenMode{
    PostFSR = 0, DressedDRp1, DressedDRp4, PreFSR
};

enum GenPID{
    TOP = 6, ELECTRON = 11, MUON = 13, TAU = 15, PHOTON = 22, PROTON = 2212
};




#endif

