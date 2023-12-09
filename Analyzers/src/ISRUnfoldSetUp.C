#include "ISRUnfoldSetUp.h"

string ISRUnfoldSetUp::get_bin_name(const UnfoldSpaceName mode) {
    if (mode ==  UnfoldSpaceName::unfolded){
        return unfolded_bin->get_bin_name();
    } else {
        return folded_bin->get_bin_name();
    }
}

string ISRUnfoldSetUp::get_raw_bin_name(const UnfoldSpaceName mode) {
    if (mode ==  UnfoldSpaceName::unfolded){
        return unfolded_bin->get_raw_bin_name();
    } else {
        return folded_bin->get_raw_bin_name();

    }
}

string ISRUnfoldSetUp::get_passed_winow_name(const UnfoldSpaceName mode) {
    int passed_dimass_cut_index = -1;
    if (mode ==  UnfoldSpaceName::unfolded){
        passed_dimass_cut_index = passed_dimass_index_unfolded;
    }
    else {
        passed_dimass_cut_index = passed_dimass_index_folded;
    }
    string low_mass = to_string(dimass_cuts[passed_dimass_cut_index]);
    string high_mass = to_string(dimass_cuts[passed_dimass_cut_index+1]); 
    low_mass = low_mass.substr(0, low_mass.find('.') + 2);
    high_mass = high_mass.substr(0, high_mass.find('.') + 2);  
    return low_mass + "to" + high_mass;
}

bool ISRUnfoldSetUp::pass_reco_cuts(double reco_pt, double reco_mass){
    // for 2D
    if (is_2d){
        if (reco_pt >= dipt_low_cut && reco_pt < dipt_high_cut &&
                reco_mass >= dimass_low_cut && reco_mass < dimass_high_cut) {
            return true;
        } else {
            return false;
        }
    }
    else {
        for (unsigned int i = 0; i < dimass_cuts.size()-1; i++){
            if (reco_pt >= dipt_low_cut && reco_pt < dipt_high_cut &&
                reco_mass >= dimass_cuts[i] && reco_mass < dimass_cuts[i+1]) {
                passed_dimass_index_folded = i;
                return true;
            }
        }
        passed_dimass_index_folded = -1;
        // TODO reco_dimass_cut_passed 
        // reco_dipt_cut_passed
        // passed_dimass_index_folded = i
        return false;  // may be better to save as a property
    }
}

bool ISRUnfoldSetUp::pass_gen_cuts(double gen_pt, double gen_mass){
    if (is_2d) {
        if (gen_pt >= dipt_low_cut && gen_pt < dipt_high_cut &&
                gen_mass >= dimass_low_cut && gen_mass < dimass_high_cut) {
            return true;
        } else {
            return false;
        }
    }
    else {
        for (unsigned int i = 0; i < dimass_cuts.size()-1; i++){
            if (gen_pt >= dipt_low_cut && gen_pt < dipt_high_cut &&
                    gen_mass >= dimass_cuts[i] && gen_mass < dimass_cuts[i+1]) {
                passed_dimass_index_unfolded = i;
                return true;
            }
        }
        passed_dimass_index_unfolded = -1; 
        return false;
    }
}

bool ISRUnfoldSetUp::is_fake(double reco_pt, double reco_mass, 
       double gen_pt, double gen_mass) {

    bool reco_passed = pass_reco_cuts(reco_pt, reco_mass);
    if (reco_passed && (passed_dimass_index_folded != passed_dimass_index_unfolded))
        return true;
    else
        return false;
} 

bool ISRUnfoldSetUp::pass_reco_gen_cuts(double reco_pt, double reco_mass,
        double gen_pt, double gen_mass) {

    bool reco_passed = pass_reco_cuts(reco_pt, reco_mass);
    bool gen_passed = pass_gen_cuts(gen_pt, gen_mass);
    if (passed_dimass_index_folded == passed_dimass_index_unfolded && reco_passed && gen_passed) 
        return true;
    else 
        return false;
}
