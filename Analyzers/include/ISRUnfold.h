#ifndef ISRUnfold_h
#define ISRUnfold_h

#include <tuple>
#include <regex>
#include "SMPAnalyzerCore.h"
#include "TUnfoldBinning.h"
#include "TVectorD.h"
#include "ISRUnfoldSetUp.h"
#include "ISRUnfoldBin.h"


class ISRUnfold : public SMPAnalyzerCore{
public:

    void initializeISRUnfold(int job_num){
        if (job_num==0) write_bins=true;
        else write_bins=false;
    }
   
    void set_base_parameter(Parameter& p, bool set_reco=true);
    void set_reco_leptons(Particle& l0, Particle& l1);
    void set_gen_leptons(Particle& l0, Particle& l1);
    void set_reco_weights(map<TString,double>& weights);
    void set_gen_weights(map<TString,double>& weights);
    void set_phase_name(const UnfoldSpaceName mode, string name);

    bool pass_lepton_cuts(const UnfoldSpaceName mode);

    string get_phase_name(const UnfoldSpaceName mode);
    double get_value(const UnfoldSpaceName mode, string var_name);

    ISRUnfoldBin* create_2d_unfold_bin(string axis1_name, string bin_name, bool axis1_uf, bool axis1_of, 
            string axis2_name, string window_name, bool axis2_uf, bool axis2_of);
    void create_2d_unfold_set(ISRUnfoldBin* bin1, ISRUnfoldBin* bin2, 
            double dipt_low_cut, double dipt_high_cut, double dimass_low_cut, double dimass_high_cut, bool turn_off_bin1=false, bool turn_off_bin2=false);
    ISRUnfoldBin* create_1d_unfold_bin(string axis_name, string bin_name, bool uf, bool of);
    void create_1d_unfold_set(ISRUnfoldBin* bin1, ISRUnfoldBin* bin2, 
            double dipt_low_cut, double dipt_high_cut, vector<double>& dimass_cuts, bool turn_off_bin1=false, bool turn_off_bin2=false);
    void create_1d_unfold_set(ISRUnfoldBin* bin1, ISRUnfoldBin* bin2, 
            double dipt_low_cut, double dipt_high_cut, 
            double dimass_low_cut, double dimass_high_cut, bool turn_off_bin1=false, bool turn_off_bin2=false);
    void create_1d_unfold_set(ISRUnfoldBin* bin1, ISRUnfoldBin* bin2, 
            double dipt_low_cut, double dipt_high_cut);

    int get_bin_index(TUnfoldBinning* bin, const UnfoldSpaceName mode) const;

    void fill_unfold_fake_hist(ISRUnfoldSetUp* unfold_setup, bool apply_lepton_cuts=true);
    void fill_unfold_acceptance_hist(ISRUnfoldSetUp* unfold_setup);
    void fill_unfold_reco_hist(ISRUnfoldSetUp* unfold_setup);
    void fill_unfold_gen_hist(ISRUnfoldSetUp* unfold_setup);

    void fill_unfold_hist(ISRUnfoldSetUp* unfold_setup, const UnfoldSpaceName mode);
    void fill_unfold_hist(ISRUnfoldSetUp* unfold_setup, TString suf, Double_t weight, const UnfoldSpaceName mode);
    //void fill_1d_hists(Parameter &p, string first_axis_bin_name, TUnfoldBinning* bin, TString suf, Double_t weight, const UnfoldSpaceName mode);

    void fill_unfold_matrixs(bool apply_lepton_cuts=true); 
    void fill_unfold_fake_hists(bool apply_lepton_cuts=true); 
    void fill_unfold_gen_hists(); 
    void fill_unfold_reco_hists(); 

    void fill_unfold_response_matrix(ISRUnfoldSetUp* unfold_setup, bool apply_lepton_cuts=true);
    void fill_unfold_response_matrix(ISRUnfoldSetUp* unfold_setup, TString suf, Double_t reco_weight, Double_t gen_weight);
    //void fill_1d_response_matrixs(Parameter &p, TString folded_bin_name, TString unfolded_bin_name, TString suf, Double_t reco_weight, Double_t gen_weight);

    virtual void WriteHist();

    ISRUnfold();
    ~ISRUnfold();
    
private :
    Parameter p;  //
    Particle reco_isr_l0, reco_isr_l1;
    Particle gen_isr_l0, gen_isr_l1;
    double reco_dipt, reco_dimass;
    double gen_dipt, gen_dimass;
    string reco_phase_name, gen_phase_name;
    map<TString,double> reco_weights;
    map<TString,double> gen_weights;

    map<TString, TUnfoldBinning*> map_unfold_2d_bins; 
    map<TString, TUnfoldBinning*> map_unfold_1d_bins; 
    bool write_bins=false;

    vector<ISRUnfoldSetUp*> unfold_setups;
};
    
#endif

