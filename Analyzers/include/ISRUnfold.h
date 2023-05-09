#ifndef ISRUnfold_h
#define ISRUnfold_h

#include <tuple>
#include "SMPAnalyzerCore.h"
#include "TUnfoldBinning.h"

enum class TUnfold_Bin{
    
    folded_bin=0, unfolded_bin, response_matrix
    
};

class ISRUnfold : public SMPAnalyzerCore{
public:

    class TUnfoldParameter{

    public:
        
        TUnfoldParameter(){}
        ~TUnfoldParameter(){}
       
        // parameter for 2D binning  
        TUnfoldParameter(const vector<double>& axis0_folded,
                         const vector<double>& axis0_unfolded,
                         const vector<double>& axis1_folded,
                         const vector<double>& axis1_unfolded,
                         bool axis0_uf, bool axis0_of, bool axis1_uf, bool axis1_of,
                         const char* axis0_var_name, const char* axis1_var_name, 
                         const char* axis0_folded_bin_name, const char* axis0_unfolded_bin_name,
                         const char* axis1_folded_bin_name, const char* axis1_unfolded_bin_name) : 

            first_axis_folded(axis0_folded), 
            first_axis_unfolded(axis0_unfolded),
            second_axis_folded(axis1_folded), 
            second_axis_unfolded(axis1_unfolded),
            use_first_axis_uf{axis0_uf}, use_first_axis_of{axis0_of}, use_second_axis_uf{axis1_uf}, use_second_axis_of{axis1_of},
            first_axis_var_name{axis0_var_name}, second_axis_var_name{axis1_var_name}, 
            first_axis_folded_bin_name{axis0_folded_bin_name}, first_axis_unfolded_bin_name{axis0_unfolded_bin_name},
            second_axis_folded_bin_name{axis1_folded_bin_name}, second_axis_unfolded_bin_name{axis1_unfolded_bin_name},
            is_2D{true} 
        {
                                                                                          n_first_axis_folded    = first_axis_folded.size()-1; // number of bins
                                                                                          n_first_axis_unfolded  = first_axis_unfolded.size()-1;
                                                                                          n_second_axis_folded   = second_axis_folded.size()-1;
                                                                                          n_second_axis_unfolded = second_axis_unfolded.size()-1;

                                                                                          var_name = "["+(string)first_axis_var_name+"-"+(string)second_axis_var_name+"]";
        }

        /*
        // parameter for 1D binning  
        TUnfoldParameter(const int n_axis0_folded, const double* axis0_folded,
                         const int n_axis0_unfolded, const double* axis0_unfolded,
                         bool first_axis_uf, bool first_axis_of,
                         const char* first_axis_var_name) : n_first_axis_folded{n_axis0_folded}, first_axis_folded{axis0_folded}, n_first_axis_unfolded{n_axis0_unfolded}, first_axis_unfolded{axis0_unfolded},
        use_first_axis_uf{first_axis_uf}, use_first_axis_of{first_axis_of}, 
        first_axis_var_name{first_axis_var_name}, is_2D{false} {
            bin_name = "1D_"+(string)first_axis_var_name;
        }
        */
        
        int n_first_axis_folded;
        int n_first_axis_unfolded;
        const vector<double> first_axis_folded;
        const vector<double> first_axis_unfolded;
        
        int n_second_axis_folded;
        int n_second_axis_unfolded;
        const vector<double> second_axis_folded;
        const vector<double> second_axis_unfolded;
        
        bool use_first_axis_uf;
        bool use_first_axis_of;
        bool use_second_axis_uf;
        bool use_second_axis_of;
        
        string var_name;
        const char* first_axis_var_name;
        const char* second_axis_var_name;

        const char* first_axis_folded_bin_name;
        const char* first_axis_unfolded_bin_name;
        const char* second_axis_folded_bin_name;
        const char* second_axis_unfolded_bin_name;
        
        bool is_2D=false;

    };
   
    void initializeISRUnfold(int job_num){
        if (job_num==0) write_bins=true;
        else write_bins=false;
    }
    
    TUnfoldBinning* get_bin_pointer(TString channelname, const TUnfoldParameter& par, const TUnfold_Bin mode);
    void fill_unfold_hists(TString channelname, TString pre, TString suf, 
            Particle* l0, Particle* l1, map<TString,double> weights, const TUnfoldParameter& par, const TUnfold_Bin mode, const TString tunfold_prefix = "");
    void fill_unfold_response_matrixs(TString channelname, TString pre, TString suf,
            Particle* l0, Particle* l1, Particle* unfolded_l0, Particle* unfolded_l1, map<TString,double> reco_weights, map<TString,double> gen_weights, const TUnfoldParameter& par, const TString tunfold_prefix=""); // response matrix,
    
    void fill_unfold_hist(TString hname, Double_t value, map<TString,double> weights, TUnfoldBinning* bin_pointer = nullptr, bool is_2D = true);
    void fill_unfold_hist(TString hname, Double_t value, Double_t weight, TUnfoldBinning* bin_pointer = nullptr, bool is_2D = true);
    void fill_unfold_response_matrix(TString hname, Double_t value_folded, Double_t value_unfolded, map<TString,double> reco_weights, map<TString,double> gen_weights, TUnfoldBinning* bin_pointer_folded=nullptr, TUnfoldBinning* bin_pointer_unfolded=nullptr, bool is_2D = true);
    void fill_unfold_response_matrix(TString hname, Double_t value_folded, Double_t value_unfolded, Double_t reco_weight, Double_t gen_weight, TUnfoldBinning* bin_pointer_folded=nullptr, TUnfoldBinning* bin_pointer_unfolded=nullptr, bool is_2D = true);
    bool is_same_bin(TString channelname, TString pre, Double_t variable1, Double_t variabl2, const TUnfoldParameter& par, const TUnfold_Bin mode);
    
    virtual void WriteHist();

    // void set_tunfold_binnings();
    
    ISRUnfold();
    ~ISRUnfold();
    
    private :
   
    map<TString, TUnfoldBinning*> map_tunfoldbins;
    bool write_bins=false;
};
    
#endif

