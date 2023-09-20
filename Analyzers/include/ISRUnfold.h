#ifndef ISRUnfold_h
#define ISRUnfold_h

#include <tuple>
#include <regex>
#include "SMPAnalyzerCore.h"
#include "TUnfoldBinning.h"
#include "TVectorD.h"

enum class TUnfoldBin{
    folded_bin=0, 
    unfolded_bin, 
    response_matrix
};

const map<TString, map<TString, vector<double>>> bins = {
    {"dipt", { 
                 {"coarse", {0., 4., 8., 12., 18., 28., 40., 55., 75., 100.}},
                 {"fine", {0., 2., 4., 6., 8., 10., 12., 14., 18., 23, 28., 34., 40., 47.5, 55., 65., 75., 87.5, 100.}},
                 {"coarse_extended", {0., 2., 4., 6., 8., 10., 12., 14., 18., 22, 28., 37., 52., 80., 130., 190., 250., 350, 450., 1000.}},
                 {"fine_extended", {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 16., 18., 20., 22., 25., 28., 32., 37., 43., 52., 65., 80., 100., 130., 160., 190., 220., 250., 300., 350., 400., 450., 500., 1000.}},
             }
    },
    {"dimass", {
                 {"coarse", {55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000}},
                 {"fine", {55, 57.5, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000}}
               }
    }
};

const map<TString, map<TString, vector<double>>> windows = {
    {"dipt", {
                 {"window_v1", {0, 100}},
                 {"window_v2", {0, 1000}}
             }
    },
    {"dimass", {
                   {"window_v1", {55., 64., 81., 101., 200., 1000.}}
               }
    }
};

class ISRUnfold : public SMPAnalyzerCore{
public:

    // explicit 2D
    // TODO get 1D from 2D definition
    // TODO effect of underflow/overflow cut
    class TUnfoldParameter{
    public:
        TUnfoldParameter(){}
        ~TUnfoldParameter(){}
       
        // parameter for 2D binning  
        TUnfoldParameter(const vector<double>& axis0_bin,
                         const vector<double>& axis1_bin,
                         bool axis0_uf, bool axis0_of, bool axis1_uf, bool axis1_of,
                         string axis0_var_name, string axis1_var_name, 
                         string axis0_bin_name, string axis1_bin_name): 

            first_axis{axis0_bin}, second_axis{axis1_bin},
            use_first_axis_uf{axis0_uf}, use_first_axis_of{axis0_of}, use_second_axis_uf{axis1_uf}, use_second_axis_of{axis1_of},
            first_axis_var_name{axis0_var_name}, second_axis_var_name{axis1_var_name}, 
            first_axis_bin_name{axis0_bin_name}, second_axis_bin_name{axis1_bin_name}
        {
            n_first_axis = first_axis.size()-1;  // number of bins
            n_second_axis = second_axis.size()-1;

            var_name = "["+first_axis_var_name+"-"+second_axis_var_name+"]";
            if(use_first_axis_uf && use_first_axis_of){ 
                first_axis_bin_name += "_UO";
            }
            else if(use_first_axis_uf){ 
                first_axis_bin_name += "_U";
            }
            else if(use_first_axis_of){
                first_axis_bin_name += "_O";
            }

            if(use_second_axis_uf && use_second_axis_of){ 
                second_axis_bin_name += "_UO";
            }
            else if(use_second_axis_uf){ 
                second_axis_bin_name += "_U";
            }
            else if(use_second_axis_of){
                second_axis_bin_name += "_O";
            }
        }

        int n_first_axis;
        int n_second_axis;
        const vector<double> first_axis;
        const vector<double> second_axis;
        
        bool use_first_axis_uf;
        bool use_first_axis_of;
        bool use_second_axis_uf;
        bool use_second_axis_of;
        
        string var_name;
        string first_axis_var_name;
        string second_axis_var_name;

        string first_axis_bin_name;
        string second_axis_bin_name;
        
        TUnfoldBinning* create_2d_tunfold_bin(const TUnfoldBin mode) const;
        string get_bin_name(const TUnfoldBin mode) const;
    };
   
    void initializeISRUnfold(int job_num){
        if (job_num==0) write_bins=true;
        else write_bins=false;
    }
   
    void set_pt_mass(const TUnfoldBin mode, double pt, double mass); 
    void set_phase_name(const TUnfoldBin mode, string name);
    string get_phase_name(const TUnfoldBin mode);
    double get_value(const TUnfoldBin mode, string var_name);

    void create_2d_folded_bin(string axis1_name, string bin_name, bool axis1_uf, bool axis1_of, 
            string axis2_name, string window_name, bool axis2_uf, bool axis2_of, bool hist_on=false);
    void create_2d_unfolded_bin(string axis1_name, string bin_name, bool axis1_uf, bool axis1_of, 
            string axis2_name, string window_name, bool axis2_uf, bool axis2_of, bool hist_on=false);

    int get_bin_index(TUnfoldBinning* bin, const TUnfoldBin mode) const;

    void fill_unfold_hists(Parameter &p, Particle* l0, Particle* l1, map<TString,double> weights, const TUnfoldBin mode, TString phase_name = "");
    void fill_unfold_hist(Parameter &p, map<TString,double> weights, const TUnfoldBin mode);
    void fill_unfold_hist(Parameter &p, TString suf, Double_t weight, const TUnfoldBin mode);
    void fill_1d_hists(Parameter &p, string first_axis_bin_name, TUnfoldBinning* bin, TString suf, Double_t weight, const TUnfoldBin mode);

    void fill_unfold_response_matrixs(Parameter &p, Particle* l0, Particle* l1, Particle* unfolded_l0, Particle* unfolded_l1,
            map<TString,double> reco_weights, map<TString,double> gen_weights, TString reco_phase_name="", TString gen_phase_name=""); // response matrix,
    void fill_unfold_response_matrix(Parameter &p, map<TString,double> reco_weights, map<TString,double> gen_weights);
    void fill_unfold_response_matrix(Parameter &p, TString suf, Double_t reco_weight, Double_t gen_weight);
    void fill_1d_response_matrixs(Parameter &p, TString folded_bin_name, TString unfolded_bin_name, TString suf, Double_t reco_weight, Double_t gen_weight);
    
    virtual void WriteHist();

    // void set_tunfold_binnings();
    ISRUnfold();
    ~ISRUnfold();
    
    private :
    double reco_pt, reco_mass;
    double gen_pt, gen_mass;
    string reco_phase_name, gen_phase_name;
    map<TString, TUnfoldBinning*> map_folded_bins;
    map<TString, TUnfoldBinning*> map_unfolded_bins; 
    map<TString, bool> map_folded_bin_flags; 
    map<TString, bool> map_unfolded_bin_flags; 
    bool write_bins=false;
};
    
#endif

