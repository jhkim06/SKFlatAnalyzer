#ifndef ISRUnfoldBin_h
#define ISRUnfoldBin_h

#include <tuple>
#include <regex>
#include "SMPAnalyzerCore.h"
#include "TUnfoldBinning.h"
#include "TVectorD.h"

const map<TString, map<TString, vector<double>>> bins = {
    {"dipt", { 
                 {"coarse", {0., 4., 8., 12., 18., 28., 40., 55., 75., 100.}},
                 {"coarse_v3", {0., 6., 12., 18., 28., 40., 55., 75., 100.}},
                 {"fine", {0., 2., 4., 6., 8., 10., 12., 14., 18., 23, 28., 34., 40., 47.5, 55., 65., 75., 87.5, 100.}},

                 //{"coarse", {0., 4., 8., 12., 17., 24.,  32.,  43.,  65.,  100}},
                 //{"coarse_v3", {0., 6., 12., 20,  32.,  43.,  65.,  100}},
                 //{"fine", {0., 2., 4., 6., 8., 10., 12., 14., 17., 20., 24., 28., 32., 37., 43., 52., 65., 80., 100}},

                 //{"coarse_extended", {0., 4., 8., 12., 17., 24.,  32.,  43.,  65.,  100.,  150.,  210.,  300.,  400., 500., 1000.}},
                 //{"coarse_v3_extended", {0., 6., 12., 18., 24.,  32.,  43.,  65.,  100.,  150.,  210.,  300.,  400., 500., 1000.}},
                 //{"fine_extended", {0., 2., 4., 6., 8., 10., 12., 14., 17., 20., 24., 28., 32., 37., 43., 52., 65., 80., 100., 125., 150., 180., 210., 250., 300., 350., 400., 450., 500., 1000.}},

                 {"coarse_extended", {0., 4., 8., 12., 17., 24.,  32.,  43.,  65.,  100.,  150.,  210.,  300.,  400., 500.}},
                 {"coarse_v3_extended", {0., 6., 12., 20,  32.,  43.,  65.,  100.,  150.,  210.,  300.,  400., 500.}},
                 {"fine_extended", {0., 2., 4., 6., 8., 10., 12., 14., 17., 20., 24., 28., 32., 37., 43., 52., 65., 80., 100., 125., 150., 180., 210., 250., 300., 350., 400., 450., 500.}},
             }
    },
    {"dimass", {
                   {"coarse", {55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000}},
                   {"coarse_high_mass_v1", {106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220}},
                   {"coarse_high_mass_v2", {110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243}},
                   {"coarse_high_mass_v3", {115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273}},
                   {"coarse_high_mass_v4", {120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320}},
                   {"coarse_high_mass_v5", {126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380}},
                   {"fine", {55, 57.5, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78.5, 81, 83.5, 86, 88.5, 91, 93.5, 96, 98.5, 101, 103.5, 106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380, 410, 440, 475, 510, 555, 600, 650, 700, 765, 830, 915, 1000}},
                   {"fine_high_mass_v1", {106, 108, 110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220}},
                   {"fine_high_mass_v2", {110, 112.5, 115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243}},
                   {"fine_high_mass_v3", {115, 117.5, 120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273}},
                   {"fine_high_mass_v4", {120, 123, 126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320}},
                   {"fine_high_mass_v5", {126, 129.5, 133, 137, 141, 145.5, 150, 155, 160, 165.5, 171, 178, 185, 192.5, 200, 210, 220, 231.5, 243, 258, 273, 296.5, 320, 350, 380}}
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

class ISRUnfoldBin{
public:

    ISRUnfoldBin(){}
    ~ISRUnfoldBin(){}
    // for 2D binning
    ISRUnfoldBin(string axis0_var_name, string axis1_var_name, 
                 string axis0_bin_name, string axis1_bin_name, 
                 bool axis0_uf, bool axis0_of, bool axis1_uf, bool axis1_of):
        first_axis_var_name{axis0_var_name}, second_axis_var_name{axis1_var_name}, 
        first_axis_original_bin_name{axis0_bin_name}, second_axis_original_bin_name{axis1_bin_name},
        use_first_axis_uf{axis0_uf}, use_first_axis_of{axis0_of}, 
        use_second_axis_uf{axis1_uf}, use_second_axis_of{axis1_of}
    {
        is_2d = true;
        var_name = "["+first_axis_var_name+"-"+second_axis_var_name+"]";

        first_axis_bin_name = first_axis_original_bin_name + add_uo_postfix(use_first_axis_uf, use_first_axis_of);
        second_axis_bin_name = second_axis_original_bin_name + add_uo_postfix(use_second_axis_uf, use_second_axis_of);
    }
    // for 1D
    ISRUnfoldBin(string axis_var_name, string axis_bin_name, bool uf=false, bool of=false):
        first_axis_var_name{axis_var_name}, first_axis_original_bin_name{axis_bin_name},
        use_first_axis_uf{uf}, use_first_axis_of{of}
    {
        is_2d = false;
        var_name = first_axis_var_name;

        first_axis_bin_name = first_axis_original_bin_name + add_uo_postfix(use_first_axis_uf, use_first_axis_of);
        if (first_axis_var_name == "dipt")
            second_axis_var_name = "dimass";
        else 
            second_axis_var_name = "dipt";
    }
    TUnfoldBinning* create_2d_tunfold_bin() const;
    TUnfoldBinning* create_1d_bin();
    vector<double> get_boundary_bin_edges();
    string get_bin_name() const;
    inline string get_var_name() {
        return var_name;
    }
    inline string get_first_axis_var_name() {
        return first_axis_var_name;
    }
    inline string get_second_axis_var_name() {
        return second_axis_var_name;
    }
    string get_raw_bin_name();
private:

    bool is_2d;
    string add_uo_postfix(bool use_axis_uf, bool use_axis_of);
    string first_axis_var_name;
    string second_axis_var_name;
    string first_axis_original_bin_name;
    string second_axis_original_bin_name;
    string first_axis_bin_name;
    string second_axis_bin_name;

    bool use_first_axis_uf;
    bool use_first_axis_of;
    bool use_second_axis_uf;
    bool use_second_axis_of;

    string var_name;
};

#endif

