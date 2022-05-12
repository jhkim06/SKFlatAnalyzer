#ifndef ISRUnfold_h
#define ISRUnfold_h

#include "SMPAnalyzerCore.h"
#include "TUnfoldBinning.h"

// separate header
class ISRUnfold : public SMPAnalyzerCore{ // ISRAnalyzerHelper :: SMPAnalyzerCore
public:
    
    //TH1D* get_hist1D(TString histname);
    // write histograms
    
    // fill
    // check axis0 name
    // if 2D used , get bin index
    void fill_unfold_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> weights, bool use_coarse_hist);
    void fill_unfold_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, Particle* truth_l0, Particle* truth_l1, map<TString,double> weights);
    
    void fill_unfold_hist(TString histname, Double_t value, map<TString,double> weights);
    void fill_unfold_hist(TString histname, Double_t smeared_value, Double_t truth_value, Double_t weight, Double_t gen_weight); // response matrix, for bin zero gen*(1-rec)
    
    // for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value,weight,n_bin,x_min,x_max);
    
    void fill_unfold_hist(TString hname, Double_t value, Double_t weight);
    
    // 1D
   
    // void fill_unfold_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> weights, TUnfold_Parameters);
    // struct TUnfold_Prameters
    // void fill_tunfold_hists
    
    // assume the same bin definition used for axis 2 (so, only second_axis used)
    void create_tunfold_hist(const int n_first_axis_fine, const double* first_axis_fine,
                             const int n_first_axis_coarse, const double* first_axis_coarse,
                             const int n_second_axis, const double* second_axis,
              bool use_axis0_uf, bool use_axis0_of, bool use_axis1_uf, bool use_axis1_of,
              const char* axis0_name, const char* axis1_name) {
        
        if(!first_axis_fine || !first_axis_coarse || !second_axis){
            cout << "Check bin definitions for unfolding.\n";
            exit(EXIT_FAILURE);
        }
        is_2D = true;
        
        coarse_bin = new TUnfoldBinning("coarse");
        fine_bin = new TUnfoldBinning("fine");
        
        coarse_bin->AddAxis(axis0_name, n_first_axis_coarse, first_axis_coarse, use_axis0_uf, use_axis0_of);
        coarse_bin->AddAxis(axis1_name, n_second_axis, second_axis, use_axis1_uf, use_axis1_of);
        
        fine_bin->AddAxis(axis0_name, n_first_axis_fine, first_axis_fine, use_axis0_uf, use_axis0_of);
        fine_bin->AddAxis(axis1_name, n_second_axis, second_axis, use_axis1_uf, use_axis1_of);
    }

    ISRUnfold();
    ~ISRUnfold();
    
    private :
   
    // parameters
    bool is_2D = false;
    TUnfoldBinning* coarse_bin = nullptr;
    TUnfoldBinning* fine_bin = nullptr;
    
    //std::map< TString, TH1D* > maphistunfold_TH1D;
    //std::map< TString, TH2D* > maphistunfold_TH2D;
};
    
#endif

