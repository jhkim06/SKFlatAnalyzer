#include "ISRUnfold.h"

ISRUnfold::ISRUnfold(){
}

ISRUnfold::~ISRUnfold(){
}

vector<string> get_reg_matchs(string source_string, string re_string="\\[([^\\]]+)\\]"){
    //std::regex re("\\[([^\\]]+)\\]");
    std::regex re(re_string);
    std::smatch match;
    vector<string> matchs;
    while (regex_search(source_string, match, re)) {
        matchs.push_back(match[1].str()); 
        // suffix to find the rest of the string.
        source_string = match.suffix().str();
    }
    return matchs;
}

void ISRUnfold::set_pt_mass(const TUnfoldBin mode, double pt, double mass) {
    if (mode == TUnfoldBin::unfolded_bin){
        this->gen_pt = pt;
        this->gen_mass = mass;
    } else {
        this->reco_pt = pt;
        this->reco_mass = mass;
    }
}

void ISRUnfold::set_phase_name(const TUnfoldBin mode, string name) {
    if (mode == TUnfoldBin::unfolded_bin){
        this->gen_phase_name = name;
    } else {
        this->reco_phase_name = name;
    }
}

string ISRUnfold::get_phase_name(const TUnfoldBin mode) {
    if (mode == TUnfoldBin::unfolded_bin){
        return this->gen_phase_name;
    } else {
        return this->reco_phase_name;
    }
}

vector<double> convert_to_vector(const TVectorD* edges){
    vector<double> edges_vector;
    for (int i = 0; i < edges->GetNoElements(); ++i) {
        edges_vector.push_back((*edges)[i]);
    }
    return edges_vector;
}

string ISRUnfold::TUnfoldParameter::get_bin_name(const TUnfoldBin mode) const{

    string bin_name;
    bin_name = this->var_name + "_[" + this->first_axis_bin_name + "-" + this->second_axis_bin_name + "]";
    return bin_name;
}

TUnfoldBinning* ISRUnfold::TUnfoldParameter::create_2d_tunfold_bin(const TUnfoldBin mode) const {

    string bin_name = this->get_bin_name(mode);  // [dipt-dimass]_[fine-window_v1]
    TUnfoldBinning* bin = new TUnfoldBinning((bin_name).c_str());
    bin->AddAxis(
            this->first_axis_var_name.c_str(),
            this->n_first_axis,
            this->first_axis.data(),
            this->use_first_axis_uf,
            this->use_first_axis_of);
    bin->AddAxis(
            this->second_axis_var_name.c_str(),
            this->n_second_axis,
            this->second_axis.data(),
            this->use_second_axis_uf,
            this->use_second_axis_of);
    return bin;
}

void ISRUnfold::create_2d_folded_bin(string axis1_name, string bin_name, bool axis1_uf, bool axis1_of,
        string axis2_name, string window_name, bool axis2_uf, bool axis2_of, bool hist_on) {
    // create TUnfoldParameter
    const vector<double>& axis1 = bins.at(axis1_name).at(bin_name);
    const vector<double>& axis2 = windows.at(axis2_name).at(window_name);

    TUnfoldParameter* tunfold_par = new TUnfoldParameter(axis1, axis2, axis1_uf, axis1_of, axis2_uf, axis2_of,
            axis1_name, axis2_name, bin_name, window_name);
    string full_bin_name = tunfold_par->get_bin_name(TUnfoldBin::folded_bin);
    map_folded_bins[full_bin_name] = tunfold_par->create_2d_tunfold_bin(TUnfoldBin::folded_bin);
    map_folded_bin_flags[full_bin_name] = hist_on;
}

void ISRUnfold::create_2d_unfolded_bin(string axis1_name, string bin_name, bool axis1_uf, bool axis1_of,
        string axis2_name, string window_name, bool axis2_uf, bool axis2_of, bool hist_on) {
    const vector<double>& axis1 = bins.at(axis1_name).at(bin_name);
    const vector<double>& axis2 = windows.at(axis2_name).at(window_name);

    TUnfoldParameter* tunfold_par = new TUnfoldParameter(axis1, axis2, axis1_uf, axis1_of, axis2_uf, axis2_of,
            axis1_name, axis2_name, bin_name, window_name);
    string full_bin_name = tunfold_par->get_bin_name(TUnfoldBin::unfolded_bin);
    map_unfolded_bins[full_bin_name] = tunfold_par->create_2d_tunfold_bin(TUnfoldBin::unfolded_bin);
    map_unfolded_bin_flags[full_bin_name] = hist_on;
}

int ISRUnfold::get_bin_index(TUnfoldBinning* bin, const TUnfoldBin mode) const{

    double pt, mass;
    if (mode == TUnfoldBin::unfolded_bin){
        pt = this->gen_pt;
        mass = this->gen_mass;
    } else {
        pt = this->reco_pt;
        mass = this->reco_mass;
    }

    int index = -999;
    if (bin->GetDistributionAxisLabel(0)=="dipt"){
        index=bin->GetGlobalBinNumber(pt, mass);
    }
    else {
        index=bin->GetGlobalBinNumber(mass, pt);
    }
    return index;
}

double ISRUnfold::get_value(const TUnfoldBin mode, string var_name) {

    if (mode == TUnfoldBin::unfolded_bin){
        if (var_name == "dipt") {
            return this->gen_pt;
        } else if (var_name == "dimass") {
            return this->gen_mass;
        }
    } else {
        if (var_name == "dipt") {
            return this->reco_pt;
        } else if (var_name == "dimass") {
            return this->reco_mass;
        }
    }
}

void ISRUnfold::fill_unfold_response_matrixs(Parameter &p, Particle* l0, Particle* l1, Particle* unfolded_l0, Particle* unfolded_l1, 
        map<TString,double> reco_weights, map<TString,double> gen_weights, TString reco_phase_name, TString gen_phase_name)
{

    TLorentzVector dilepton_folded=(*l0)+(*l1);
    set_pt_mass(TUnfoldBin::folded_bin, dilepton_folded.Pt(), dilepton_folded.M()); 
    set_phase_name(TUnfoldBin::folded_bin, string(reco_phase_name));

    TLorentzVector dilepton_unfolded=(*unfolded_l0)+(*unfolded_l1);
    set_pt_mass(TUnfoldBin::unfolded_bin, dilepton_unfolded.Pt(), dilepton_unfolded.M()); 
    set_phase_name(TUnfoldBin::unfolded_bin, string(gen_phase_name));

    fill_unfold_response_matrix(p, reco_weights, gen_weights);
}

void ISRUnfold::fill_unfold_response_matrix(Parameter &p, map<TString,double> reco_weights, map<TString,double> gen_weights)
{
    // TODO check case where only gen_weights vary
    for (const auto& [suffix,reco_weight]:reco_weights){
        fill_unfold_response_matrix(p, suffix, reco_weight, gen_weights[suffix]);
    }
}

void ISRUnfold::fill_unfold_response_matrix(Parameter &p, TString suf, Double_t reco_weight, Double_t gen_weight)
{
    // loop over map_folded_bin and map_unfolded_bin 
    for (const auto& [unfolded_bin_name, unfolded_bin]: map_unfolded_bins) {
        for (const auto& [folded_bin_name, folded_bin]: map_folded_bins) {
            int unfolded_index = get_bin_index(unfolded_bin, TUnfoldBin::unfolded_bin);
            int folded_index = get_bin_index(folded_bin, TUnfoldBin::folded_bin);
            vector<string> unfolded_matchs = get_reg_matchs(string(unfolded_bin_name));
            vector<string> folded_matchs = get_reg_matchs(string(folded_bin_name));

            if (unfolded_matchs.at(0) != folded_matchs.at(0)) continue; 

            string hname = string(p.prefix) + string(p.hprefix) + "[tunfold-matrix]_" + "[" + unfolded_matchs.at(0) + "]_[" +
                this->reco_phase_name + "__" + folded_matchs.at(1) + "]_[" + 
                this->gen_phase_name + "__" + unfolded_matchs.at(1) + "]" + string(p.suffix) + string(suf);
            TH2D *this_hist = GetHist2D(hname);
            if (!this_hist){
                this_hist = (TH2D*) TUnfoldBinning::CreateHistogramOfMigrations(unfolded_bin, folded_bin, hname.data());
                this_hist->SetDirectory(NULL); 
                maphist_TH2D[hname] = this_hist;
            }
            this_hist->Fill(unfolded_index, folded_index, reco_weight);
            this_hist->Fill(unfolded_index, 0., gen_weight-reco_weight); // bin zero for 2D 

            // create matrix for 1D hists
            // option: wheter to use TUnfoldBinning 
        }
    }
}

void ISRUnfold::fill_unfold_hists(Parameter &p, Particle* l0, Particle* l1, map<TString,double> weights,
        const TUnfoldBin mode, TString phase_name){

    TLorentzVector dilepton = (*l0) + (*l1);

    set_pt_mass(mode, dilepton.Pt(), dilepton.M());
    set_phase_name(mode, string(phase_name));
    fill_unfold_hist(p, weights, mode);
}

void ISRUnfold::fill_unfold_hist(Parameter &p, map<TString,double> weights, const TUnfoldBin mode){
    // loop over weights
    for (const auto& [suffix,weight]:weights){
        fill_unfold_hist(p, suffix, weight, mode);
    }
}

 void ISRUnfold::fill_unfold_hist(Parameter &p, TString suf, Double_t weight, const TUnfoldBin mode){
    std::map<TString, TUnfoldBinning*>* map_bin;
    if (mode == TUnfoldBin::unfolded_bin) {
        map_bin = &map_unfolded_bins;
    }
    else {
        map_bin = &map_folded_bins;
    }
    // loop over bin map and fill
    for (const auto& [bin_name, bin]: *map_bin){
        // bin_name: [dipt-dimass]_[unfolded_fine_O-window_v1_UO]
        // "[" + dipt-dimass + "]_[" + level + "__unfolded_fine_O-window_v1_UO + "]" 
        vector<string> matchs = get_reg_matchs(string(bin_name));
        int index = get_bin_index(bin, mode);
        string phase_name = get_phase_name(mode);
        string hname = string(p.prefix) + string(p.hprefix) + "[tunfold-hist]_" + "[" + matchs.at(0) + "]_[" + 
            phase_name + "__" + matchs.at(1) + "]" + string(p.suffix) + string(suf);
        TH1D *this_hist = GetHist1D(hname);
        if (!this_hist){
            this_hist = (TH1D*) bin->CreateHistogram(hname.data());
            this_hist->SetDirectory(NULL);
            maphist_TH1D[hname] = this_hist;
        }
        this_hist->Fill(index, weight);
        fill_1d_hists(p, matchs.at(1), bin, suf, weight, mode); 
    }
}

void ISRUnfold::fill_1d_hists(Parameter &p, string first_axis_bin_name, TUnfoldBinning* bin, TString suf, Double_t weight, const TUnfoldBin mode) {
    
    vector<double> first_axis_edges = convert_to_vector(bin->GetDistributionBinning(0));
    string first_axis_var = string(bin->GetDistributionAxisLabel(0));
    vector<double> second_axis_edges = convert_to_vector(bin->GetDistributionBinning(1));
    string second_axis_var = string(bin->GetDistributionAxisLabel(1));

    string phase_name = get_phase_name(mode);
    for (unsigned int i = 0; i < second_axis_edges.size()-1; i++){
        string low_mass = to_string(second_axis_edges.at(i));
        string high_mass = to_string(second_axis_edges.at(i+1));
        low_mass = low_mass.substr(0, low_mass.find('.') + 1);
        high_mass = high_mass.substr(0, high_mass.find('.') + 1);
        string hname = string(p.prefix) + string(p.hprefix) + first_axis_var + "_[" + phase_name + "__" + first_axis_bin_name + "]_" + second_axis_var + "_" + low_mass + "to" + high_mass + 
            string(p.suffix) + string(suf);
        double first_axis_value = get_value(mode, first_axis_var);
        double second_axis_value = get_value(mode, second_axis_var);
        if (second_axis_value >= second_axis_edges.at(i) && second_axis_value < second_axis_edges.at(i+1))
            FillHist(hname, first_axis_value, weight, first_axis_edges.size()-1, first_axis_edges.data());
    }
}

void ISRUnfold::WriteHist(){

    SMPAnalyzerCore::WriteHist();
    outfile->cd();
    // loop over bin definition
    if (write_bins==true){
        for (std::map<TString,TUnfoldBinning*>::iterator mapit = map_folded_bins.begin(); mapit!=map_folded_bins.end(); mapit++){
            TString this_fullname=mapit->first;
            TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
            TString this_suffix=this_fullname(0,this_fullname.Last('/'));
            TDirectory *dir = outfile->GetDirectory(this_suffix);
            if(!dir){
                outfile->mkdir(this_suffix);
            }
            outfile->cd(this_suffix);
            mapit->second->Write("[tunfold-bin]_"+ this_name);
            outfile->cd();
        }
        for (std::map<TString,TUnfoldBinning*>::iterator mapit = map_unfolded_bins.begin(); mapit!=map_unfolded_bins.end(); mapit++){
            TString this_fullname=mapit->first;
            TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
            TString this_suffix=this_fullname(0,this_fullname.Last('/'));
            TDirectory *dir = outfile->GetDirectory(this_suffix);
            if(!dir){
                outfile->mkdir(this_suffix);
            }
            outfile->cd(this_suffix);
            mapit->second->Write("[tunfold-bin]_"+ this_name);
            outfile->cd();
        }
    }
}
