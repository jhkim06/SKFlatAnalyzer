#ifndef ISRUnfold_h
#define ISRUnfold_h

#include <tuple>
#include "SMPAnalyzerCore.h"
#include "TUnfoldBinning.h"

class ISRUnfold : public SMPAnalyzerCore{
public:

    class TUnfoldParameter{

    public:
        
        TUnfoldParameter(){}
        ~TUnfoldParameter(){}
        
        
        TUnfoldParameter(const int n_axis0_smeared, const double* axis0_smeared,
                         const int n_axis0_truth, const double* axis0_truth,
                         const int n_axis1_smeared, const double* axis1_smeared,
                         const int n_axis1_truth, const double* axis1_truth,
                         bool first_axis_uf, bool first_axis_of, bool second_axis_uf, bool second_axis_of,
                         const char* bin_name, const char* first_axis_name, const char* second_axis_name) : n_first_axis_smeared{n_axis0_smeared}, first_axis_smeared{axis0_smeared}, n_first_axis_truth{n_axis0_truth}, first_axis_truth{axis0_truth},
        n_second_axis_smeared{n_axis1_smeared}, second_axis_smeared{axis1_smeared}, n_second_axis_truth{n_axis1_truth}, second_axis_truth{axis1_truth},
        use_first_axis_uf{first_axis_uf}, use_first_axis_of{first_axis_of}, use_second_axis_uf{second_axis_uf}, use_second_axis_of{second_axis_of},
        bin_name{bin_name}, first_axis_name{first_axis_name}, second_axis_name{second_axis_name}, is_2D{true} { }
        
        
        int n_first_axis_smeared;
        const double* first_axis_smeared;
        int n_first_axis_truth;
        const double* first_axis_truth;
        
        int n_second_axis_smeared;
        const double* second_axis_smeared;
        int n_second_axis_truth;
        const double* second_axis_truth;
        
        bool use_first_axis_uf;
        bool use_first_axis_of;
        bool use_second_axis_uf;
        bool use_second_axis_of;
        
        const char* bin_name;
        const char* first_axis_name;
        const char* second_axis_name;
        
        bool is_2D = false;
    };
    
    //TH1D* get_hist1D(TString histname);
    // write histograms
    
    // fill
    // check axis0 name
    // if 2D used, get bin index
    void fill_unfold_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> weights, const TUnfoldParameter& par, const int mode);
    void fill_unfold_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, Particle* truth_l0, Particle* truth_l1, map<TString,double> weights);
    
    void fill_unfold_hist(TString histname, Double_t value, map<TString,double> weights, TUnfoldBinning* bin_pointer = nullptr);
    void fill_unfold_hist(TString histname, Double_t smeared_value, Double_t truth_value, Double_t weight, Double_t gen_weight); // response matrix, for bin zero gen*(1-rec)
    
    // for(const auto& [suffix,weight]:weights) FillHist(histname+suffix,value,weight,n_bin,x_min,x_max);
    
    void fill_unfold_hist(TString hname, Double_t value, Double_t weight, TUnfoldBinning* bin_pointer = nullptr);
    
    // 1D
   
    // struct TUnfold_Prameters
    // void fill_tunfold_hists
    
    // assume the same bin definition used for axis 2 (so, only second_axis used)
        
    ISRUnfold();
    ~ISRUnfold();
    
    private :
   
    map<TString, tuple<TUnfoldBinning*, TUnfoldBinning*>> map_tunfoldbins;
};
    
#endif

