#ifndef ISRUnfold_h
#define ISRUnfold_h

#include <tuple>
#include "SMPAnalyzerCore.h"
#include "TUnfoldBinning.h"

enum class TUnfold_Bin{
    
    smeared_bin=0, truth_bin, response_matrix
    
};

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
                         const char* first_axis_name, const char* second_axis_name) : n_first_axis_smeared{n_axis0_smeared}, first_axis_smeared{axis0_smeared}, n_first_axis_truth{n_axis0_truth}, first_axis_truth{axis0_truth},
        n_second_axis_smeared{n_axis1_smeared}, second_axis_smeared{axis1_smeared}, n_second_axis_truth{n_axis1_truth}, second_axis_truth{axis1_truth},
        use_first_axis_uf{first_axis_uf}, use_first_axis_of{first_axis_of}, use_second_axis_uf{second_axis_uf}, use_second_axis_of{second_axis_of},
        first_axis_name{first_axis_name}, second_axis_name{second_axis_name}, is_2D{true} {
            bin_name = "2D_"+(string)first_axis_name+"_"+(string)second_axis_name;
        }
        
        
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
        
        string bin_name;
        const char* first_axis_name;
        const char* second_axis_name;
        
        bool is_2D=false;
    };
   
    void initializeISRUnfold(int job_num){
        if (job_num==0) write_bins=true;
        else write_bins=false;
    }
    
    void fill_unfold_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> weights, const TUnfoldParameter& par, const TUnfold_Bin mode);
    void fill_unfold_response_matrixs(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, Particle* truth_l0, Particle* truth_l1, map<TString,double> reco_weights, map<TString,double> gen_weights,
                          const TUnfoldParameter& par); // response matrix,
    
    void fill_unfold_hist(TString hname, Double_t value, map<TString,double> weights, TUnfoldBinning* bin_pointer = nullptr);
    void fill_unfold_hist(TString hname, Double_t value, Double_t weight, TUnfoldBinning* bin_pointer = nullptr);
    void fill_unfold_response_matrix(TString hname, Int_t value_smeared, Int_t value_truth, map<TString,double> reco_weights, map<TString,double> gen_weights, TUnfoldBinning* bin_pointer_smeared=nullptr, TUnfoldBinning* bin_pointer_truth=nullptr);
    void fill_unfold_response_matrix(TString hname, Int_t value_smeared, Int_t value_truth, Double_t reco_weight, Double_t gen_weight, TUnfoldBinning* bin_pointer_smeared=nullptr, TUnfoldBinning* bin_pointer_truth=nullptr);
    
    virtual void WriteHist();
   
    
    ISRUnfold();
    ~ISRUnfold();
    
    private :
   
    map<TString, tuple<TUnfoldBinning*, TUnfoldBinning*>> map_tunfoldbins;
    bool write_bins=false;
};
    
#endif

