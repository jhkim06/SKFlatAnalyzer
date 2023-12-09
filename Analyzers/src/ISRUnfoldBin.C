#include "ISRUnfoldBin.h"

string ISRUnfoldBin::get_bin_name() const{

    string bin_name;
    if (is_2d)
        bin_name = this->var_name + "_[" + this->first_axis_bin_name + "-" + this->second_axis_bin_name + "]";
    else 
        bin_name = this->var_name + "_[" + this->first_axis_bin_name + "]";
    return bin_name;
}

string ISRUnfoldBin::get_raw_bin_name() {
    string bin_name;
    if (is_2d)
        bin_name = this->first_axis_bin_name + "-" + this->second_axis_bin_name;
    else 
        bin_name = this->first_axis_bin_name;
    return bin_name;
}

string ISRUnfoldBin::add_uo_postfix(bool use_axis_uf, bool use_axis_of) {  

    string postfix;
    if(use_axis_uf && use_axis_of){
        postfix = "_UO"; 
    }
    else if(use_axis_uf){ 
        postfix = "_U"; 
    }
    else if(use_axis_of){  
        postfix = "_O";
    }
    else
        postfix = "";   
    return postfix;
}

TUnfoldBinning* ISRUnfoldBin::create_2d_tunfold_bin() const {

    const vector<double>& axis1 = bins.at(first_axis_var_name).at(first_axis_original_bin_name);  
    const vector<double>& axis2 = windows.at(second_axis_var_name).at(second_axis_original_bin_name);

    int n_first_axis = axis1.size()-1;
    int n_second_axis = axis2.size()-1;

    string bin_name = this->get_bin_name();  // [dipt-dimass]_[fine-window_v1]
    TUnfoldBinning* bin = new TUnfoldBinning((bin_name).c_str()); // can get this name later??
    bin->AddAxis(
            this->first_axis_var_name.c_str(),
            n_first_axis,
            axis1.data(),
            this->use_first_axis_uf,
            this->use_first_axis_of);
    bin->AddAxis(
            this->second_axis_var_name.c_str(),
            n_second_axis,
            axis2.data(),
            this->use_second_axis_uf,
            this->use_second_axis_of);
    return bin;
}

vector<double> ISRUnfoldBin::create_1d_bin() {
    vector<double> bin = bins.at(first_axis_var_name).at(first_axis_original_bin_name);
    return bin;
}

vector<double> ISRUnfoldBin::get_boundary_bin_edges() {
    double first_edge = bins.at(first_axis_var_name).at(first_axis_original_bin_name).at(0);
    double last_edge = bins.at(first_axis_var_name).at(first_axis_original_bin_name).back();
    vector<double> edges = {first_edge, last_edge};
    return edges;
}
