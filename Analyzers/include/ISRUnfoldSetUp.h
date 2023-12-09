#ifndef ISRUnfoldSetUp_h
#define ISRUnfoldSetUp_h

#include "ISRUnfoldBin.h"

enum class UnfoldSpaceName{
    folded=0,
    unfolded
};

class ISRUnfoldSetUp {
    public:
        ISRUnfoldSetUp(){}
        ~ISRUnfoldSetUp(){}

        // 2D setup
        ISRUnfoldSetUp(ISRUnfoldBin *folded_bin, ISRUnfoldBin *unfolded_bin, 
                double dipt_low_cut, double dipt_high_cut,
                double dimass_low_cut, double dimass_high_cut, bool turn_off_bin1=false, bool turn_off_bin2=false): 

            folded_bin{folded_bin}, unfolded_bin{unfolded_bin}, 
            dipt_low_cut{dipt_low_cut}, dipt_high_cut{dipt_high_cut}, 
            dimass_low_cut{dimass_low_cut}, dimass_high_cut{dimass_high_cut}
        {
                is_2d = true;
                folded_bin_turned_off = turn_off_bin1;
                unfolded_bin_turned_off = turn_off_bin2;
        }
        // 1D
        ISRUnfoldSetUp(ISRUnfoldBin *folded_bin, ISRUnfoldBin *unfolded_bin,
                double dipt_low_cut, double dipt_high_cut,
                vector<double>& dimass_cuts_, bool turn_off_bin1=false, bool turn_off_bin2=false):

            folded_bin{folded_bin}, unfolded_bin{unfolded_bin}, 
            dipt_low_cut{dipt_low_cut}, dipt_high_cut{dipt_high_cut}
        {
            is_2d = false;
            folded_bin_turned_off = turn_off_bin1;
            unfolded_bin_turned_off = turn_off_bin2;
            dimass_cuts = dimass_cuts_;
            passed_dimass_index_folded = -1;
            passed_dimass_index_unfolded = -1;
        }
        inline bool bin_turned_off(const UnfoldSpaceName mode){
            if (mode ==  UnfoldSpaceName::unfolded) {
                return unfolded_bin_turned_off;
            } else {
                return folded_bin_turned_off;
            }
        }
        bool pass_reco_cuts(double reco_pt, double reco_mass);
        bool pass_gen_cuts(double gen_pt, double gen_mass);
        bool pass_reco_gen_cuts(double reco_pt, double reco_mass, double gen_pt, double gen_mass);
        bool is_fake(double reco_pt, double reco_mass, double gen_pt, double gen_mass);
        inline bool is_2d_unfold(){
            return is_2d;
        }

        inline string get_var_name() {return folded_bin->get_var_name();}
        inline string get_first_axis_var_name() {return folded_bin->get_first_axis_var_name();}
        inline string get_second_axis_var_name() {return folded_bin->get_second_axis_var_name();}
        string get_bin_name(const UnfoldSpaceName mode);
        string get_raw_bin_name(const UnfoldSpaceName mode);
        string get_passed_winow_name(const UnfoldSpaceName mode);
        inline string get_dipt_range() {
            string low_pt = to_string(dipt_low_cut);
            string high_pt = to_string(dipt_high_cut);
            return low_pt.substr(0, low_pt.find('.') + 2) + "to" + high_pt.substr(0, high_pt.find('.') + 2);
        } 
        // 1D setup
    private:
        ISRUnfoldBin* folded_bin;
        ISRUnfoldBin* unfolded_bin;
        bool is_2d;
        double dipt_low_cut;  // overflow/underflow cuts
        double dipt_high_cut;
        double dimass_low_cut;  // used with 2D
        double dimass_high_cut;  // used with 2D
        bool folded_bin_turned_off;
        bool unfolded_bin_turned_off;


        vector<double> dimass_cuts;
        int passed_dimass_index_folded = -1;
        int passed_dimass_index_unfolded = -1;
};
#endif
