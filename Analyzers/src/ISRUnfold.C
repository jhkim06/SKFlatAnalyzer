#include "ISRUnfold.h"

ISRUnfold::ISRUnfold(){
}

ISRUnfold::~ISRUnfold(){
}

vector<string> get_reg_matchs(string source_string, string re_string="\\[([^\\]]+)\\]"){
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

void ISRUnfold::set_base_parameter(Parameter& p, bool set_reco) {
    this->p = p;
    if (set_reco) {
        set_reco_leptons(*p.lepton0, *p.lepton1);
        set_reco_weights(p.weightmap);
    }
}

void ISRUnfold::set_reco_weights(map<TString,double>& weights) {
    reco_weights = weights;
}

void ISRUnfold::set_gen_weights(map<TString,double>& weights) {
    gen_weights = weights;
}

void ISRUnfold::set_reco_leptons(Particle& l0, Particle& l1) {
    reco_isr_l0 = l0;
    reco_isr_l1 = l1;

    reco_dipt = (l0 + l1).Pt();
    reco_dimass = (l0 + l1).M();
}

void ISRUnfold::set_gen_leptons(Particle& l0, Particle& l1) {
    gen_isr_l0 = l0;
    gen_isr_l1 = l1;

    gen_dipt = (l0 + l1).Pt();
    gen_dimass = (l0 + l1).M();
}

void ISRUnfold::set_phase_name(const UnfoldSpaceName mode, string name) {
    if (mode == UnfoldSpaceName::unfolded){
        this->gen_phase_name = name;
    } else {
        this->reco_phase_name = name;
    }
}

string ISRUnfold::get_phase_name(const UnfoldSpaceName mode) {
    if (mode == UnfoldSpaceName::unfolded){
        return this->gen_phase_name;
    } else {
        return this->reco_phase_name;
    }
}

bool ISRUnfold::pass_lepton_cuts(const UnfoldSpaceName mode) {

    Particle* lepton0;
    Particle* lepton1;
    if (mode == UnfoldSpaceName::unfolded) {
        lepton0 = &gen_isr_l0;
        lepton1 = &gen_isr_l1;
    } else {
        lepton0 = &reco_isr_l0;
        lepton1 = &reco_isr_l1;
    }

    bool pass_lepton_pt_cut = ((*lepton0).Pt() > p.c.lepton0pt && (*lepton1).Pt() > p.c.lepton1pt) || 
        ((*lepton1).Pt() > p.c.lepton0pt && (*lepton0).Pt() > p.c.lepton1pt);

    double eta_cut = p.c.lepton_max_eta;
    bool pass_lepton_eta_cut = (fabs((*lepton0).Eta()) < eta_cut) && (fabs((*lepton1).Eta()) < eta_cut);

    if (pass_lepton_pt_cut && pass_lepton_eta_cut) {
        return true;
    } else {
        return false;
    }
}

vector<double> convert_to_vector(const TVectorD* edges){
    vector<double> edges_vector;
    for (int i = 0; i < edges->GetNoElements(); ++i) {
        edges_vector.push_back((*edges)[i]);
    }
    return edges_vector;
}

ISRUnfoldBin* ISRUnfold::create_2d_unfold_bin(string axis1_name, string bin_name, bool axis1_uf, bool axis1_of,
        string axis2_name, string window_name, bool axis2_uf, bool axis2_of) {

    ISRUnfoldBin* tunfold_par = new ISRUnfoldBin(axis1_name, axis2_name, bin_name, window_name,
            axis1_uf, axis1_of, axis2_uf, axis2_of);

    string full_bin_name = tunfold_par->get_bin_name();
    map_unfold_2d_bins[full_bin_name] = tunfold_par->create_2d_tunfold_bin();

    // return parameter and use it for unfold setup
    return tunfold_par;
}

ISRUnfoldBin* ISRUnfold::create_1d_unfold_bin(string axis_name, string bin_name) {

    ISRUnfoldBin* tunfold_par = new ISRUnfoldBin(axis_name, bin_name);

    string full_bin_name = tunfold_par->get_bin_name();
    map_unfold_1d_bins[full_bin_name] = tunfold_par->create_1d_bin();

    return tunfold_par;
}

void ISRUnfold::create_2d_unfold_set(ISRUnfoldBin* bin1, ISRUnfoldBin* bin2,
        double dipt_low_cut, double dipt_high_cut, 
        double dimass_low_cut, double dimass_high_cut,
        bool turn_off_bin1, bool turn_off_bin2) {

    // parameter
    ISRUnfoldSetUp* unfold_setup = new ISRUnfoldSetUp(bin1, bin2, 
            dipt_low_cut, dipt_high_cut, 
            dimass_low_cut, dimass_high_cut,
            turn_off_bin1, turn_off_bin2);
    // add to vector
    unfold_setups.push_back(unfold_setup);
}

void ISRUnfold::create_1d_unfold_set(ISRUnfoldBin* bin1, ISRUnfoldBin* bin2,  
        double dipt_low_cut, double dipt_high_cut,
        vector<double>& dimass_cuts, bool turn_off_bin1, bool turn_off_bin2) {

    ISRUnfoldSetUp* unfold_setup = new ISRUnfoldSetUp(bin1, bin2,
            dipt_low_cut, dipt_high_cut,
            dimass_cuts, turn_off_bin1, turn_off_bin2);

    unfold_setups.push_back(unfold_setup);
}

void ISRUnfold::create_1d_unfold_set(ISRUnfoldBin* bin1, ISRUnfoldBin* bin2,  
        double dipt_low_cut, double dipt_high_cut,
        double dimass_low_cut, double dimass_high_cut, bool turn_off_bin1, bool turn_off_bin2) {

    // make vector
    vector<double> dimass_cuts = {dimass_low_cut, dimass_high_cut};
    ISRUnfoldSetUp* unfold_setup = new ISRUnfoldSetUp(bin1, bin2,
            dipt_low_cut, dipt_high_cut,
            dimass_cuts, turn_off_bin1, turn_off_bin2);

    unfold_setups.push_back(unfold_setup);
}

void ISRUnfold::create_1d_unfold_set(ISRUnfoldBin* bin1, ISRUnfoldBin* bin2,  
        double dipt_low_cut, double dipt_high_cut) {

    // make vector
    vector<double> dimass_cuts = bin1->get_boundary_bin_edges();
    ISRUnfoldSetUp* unfold_setup = new ISRUnfoldSetUp(bin1, bin2,
            dipt_low_cut, dipt_high_cut,
            dimass_cuts);

    unfold_setups.push_back(unfold_setup);
}

void ISRUnfold::fill_unfold_matrixs(bool apply_lepton_cuts) {
   
   for (const auto unfold_setup: unfold_setups){
       fill_unfold_response_matrix(unfold_setup, apply_lepton_cuts);
   }
}

void ISRUnfold::fill_unfold_fake_hists(bool apply_lepton_cuts) {
   
   for (const auto unfold_setup: unfold_setups){
       fill_unfold_fake_hist(unfold_setup, apply_lepton_cuts);
   }
}

void ISRUnfold::fill_unfold_gen_hists() {
   
   for (const auto unfold_setup: unfold_setups){
       fill_unfold_gen_hist(unfold_setup);
   }
}

void ISRUnfold::fill_unfold_reco_hists() {
   
   for (const auto unfold_setup: unfold_setups){
       fill_unfold_reco_hist(unfold_setup);
   }
}

int ISRUnfold::get_bin_index(TUnfoldBinning* bin, const UnfoldSpaceName mode) const{

    double pt, mass;
    if (mode == UnfoldSpaceName::unfolded){
        pt = this->gen_dipt;
        mass = this->gen_dimass;
    } else {
        pt = this->reco_dipt;
        mass = this->reco_dimass;
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

double ISRUnfold::get_value(const UnfoldSpaceName mode, string var_name) {

    if (mode == UnfoldSpaceName::unfolded){
        if (var_name == "dipt") {
            return this->gen_dipt;
        } else if (var_name == "dimass") {
            return this->gen_dimass;
        }
    } else {
        if (var_name == "dipt") {
            return this->reco_dipt;
        } else if (var_name == "dimass") {
            return this->reco_dimass;
        }
    }
}

void ISRUnfold::fill_unfold_response_matrix(ISRUnfoldSetUp* unfold_setup, bool apply_lepton_cuts)
{
    bool pass_reco_gen = unfold_setup->pass_reco_gen_cuts(reco_dipt, reco_dimass, gen_dipt, gen_dimass);
    bool pass_gen_lepton_cuts = pass_lepton_cuts(UnfoldSpaceName::unfolded);  // TODO ensure reco lepton cuts passed
    if (!apply_lepton_cuts) pass_gen_lepton_cuts = true;
    if (pass_reco_gen && pass_gen_lepton_cuts){ 
        // TODO check case where only gen_weights vary
        for (const auto& [suffix,reco_weight]:reco_weights){
            fill_unfold_response_matrix(unfold_setup, suffix, reco_weight, gen_weights[suffix]);
        }
    }
    else{
        return;
    }
}

void ISRUnfold::fill_unfold_response_matrix(ISRUnfoldSetUp* unfold_setup, TString suf, Double_t reco_weight, Double_t gen_weight)
{
    // get bin names
    string unfolded_bin_name = unfold_setup->get_bin_name(UnfoldSpaceName::unfolded);
    string folded_bin_name = unfold_setup->get_bin_name(UnfoldSpaceName::folded);
    bool is_2d = unfold_setup->is_2d_unfold();

    if (is_2d) {
        TUnfoldBinning* unfolded_bin = map_unfold_2d_bins[unfolded_bin_name];
        TUnfoldBinning* folded_bin = map_unfold_2d_bins[folded_bin_name];
        vector<string> unfolded_matchs = get_reg_matchs(unfolded_bin_name);
        vector<string> folded_matchs = get_reg_matchs(folded_bin_name);
        
        // get matrix_name from UnfoldSet
        string hname = string(p.prefix) + string(p.hprefix) + "[tunfold-matrix]_" + "[" + unfolded_matchs.at(0) + "]_[" +
            this->reco_phase_name + "__" + folded_matchs.at(1) + "]_[" + 
            this->gen_phase_name + "__" + unfolded_matchs.at(1) + "]" + string(p.suffix) + string(suf);
        
        TH2D *this_hist = GetHist2D(hname);
        if (!this_hist){
            this_hist = (TH2D*) TUnfoldBinning::CreateHistogramOfMigrations(unfolded_bin, folded_bin, hname.data());
            this_hist->SetDirectory(NULL); 
            maphist_TH2D[hname] = this_hist;
        }
        // get index from UnfoldSet
        int unfolded_index = get_bin_index(unfolded_bin, UnfoldSpaceName::unfolded);
        int folded_index = get_bin_index(folded_bin, UnfoldSpaceName::folded);
        this_hist->Fill(unfolded_index, folded_index, reco_weight);
        this_hist->Fill(unfolded_index, 0., gen_weight-reco_weight); // bin zero for 2D 
    }
    else {
        string var_name = unfold_setup->get_var_name();
        string second_var_range;
        if (var_name == "dipt") 
            second_var_range = unfold_setup->get_passed_winow_name(UnfoldSpaceName::folded);  // either mode will work
        else // dimass  
            second_var_range = unfold_setup->get_dipt_range();
        string hname = string(p.prefix) + string(p.hprefix) + var_name +
             "_[" + this->reco_phase_name + "__" + unfold_setup->get_raw_bin_name(UnfoldSpaceName::folded) + "]_[" +
             this->gen_phase_name + "__" + unfold_setup->get_raw_bin_name(UnfoldSpaceName::unfolded) + "]_" + 
             unfold_setup->get_second_axis_var_name() + "_" + second_var_range + string(p.suffix) + string(suf);
        
        double value_unfolded = get_value(UnfoldSpaceName::unfolded, var_name);
        double value_folded = get_value(UnfoldSpaceName::folded, var_name); 

        FillHist(hname, value_unfolded, value_folded, reco_weight, 
                map_unfold_1d_bins[unfolded_bin_name].size()-1, map_unfold_1d_bins[unfolded_bin_name].data(), 
                map_unfold_1d_bins[folded_bin_name].size()-1, map_unfold_1d_bins[folded_bin_name].data());

        // bin zero
        FillHist(hname, value_unfolded, -1, gen_weight-reco_weight, 
                map_unfold_1d_bins[unfolded_bin_name].size()-1, map_unfold_1d_bins[unfolded_bin_name].data(), 
                map_unfold_1d_bins[folded_bin_name].size()-1, map_unfold_1d_bins[folded_bin_name].data());
    }
}

void ISRUnfold::fill_unfold_fake_hist(ISRUnfoldSetUp* unfold_setup, bool apply_lepton_cuts) {
    // require passing reco and not passing gen cut
    bool pass_reco = unfold_setup->pass_reco_cuts(reco_dipt, reco_dimass);
    bool pass_fake = unfold_setup->is_fake(reco_dipt, reco_dimass, gen_dipt, gen_dimass);
    bool pass_gen_lepton_cuts = pass_lepton_cuts(UnfoldSpaceName::unfolded);
    if (!apply_lepton_cuts) pass_gen_lepton_cuts = true;

    if (pass_reco && (pass_fake || !pass_gen_lepton_cuts)){
        string original_reco_phase_name = reco_phase_name;
        reco_phase_name += "_" + gen_phase_name + "_fake";
        fill_unfold_hist(unfold_setup, UnfoldSpaceName::folded);
        reco_phase_name = original_reco_phase_name;

        // fail1: pass_fake && pass_gen_lepton_cuts
        if (pass_fake && pass_gen_lepton_cuts){
            string original_reco_phase_name = reco_phase_name;
            reco_phase_name += "_" + gen_phase_name + "_fake_fail1";
            fill_unfold_hist(unfold_setup, UnfoldSpaceName::folded);
            reco_phase_name = original_reco_phase_name;
        }
        // fail2: pass_fake && !pass_gen_lepton_cuts
        if (pass_fake && !pass_gen_lepton_cuts){
            string original_reco_phase_name = reco_phase_name;
            reco_phase_name += "_" + gen_phase_name + "_fake_fail2";
            fill_unfold_hist(unfold_setup, UnfoldSpaceName::folded);
            reco_phase_name = original_reco_phase_name;
        }
        // fail3: !pass_fake && !pass_gen_lepton_cuts
        if (!pass_fake && !pass_gen_lepton_cuts){
            string original_reco_phase_name = reco_phase_name;
            reco_phase_name += "_" + gen_phase_name + "_fake_fail3";
            fill_unfold_hist(unfold_setup, UnfoldSpaceName::folded);
            reco_phase_name = original_reco_phase_name;
        }
    }
    else {
        return;
    }
}

void ISRUnfold::fill_unfold_gen_hist(ISRUnfoldSetUp* unfold_setup) {
    bool pass_gen = unfold_setup->pass_gen_cuts(gen_dipt, gen_dimass);
    if (pass_gen){
        string original_gen_phase_name = gen_phase_name;
        gen_phase_name += "_acceptance";
        fill_unfold_hist(unfold_setup, UnfoldSpaceName::unfolded);
        gen_phase_name = original_gen_phase_name;

        bool pass_gen_lepton_cuts = pass_lepton_cuts(UnfoldSpaceName::unfolded);
        if (pass_gen_lepton_cuts){
            string original_gen_phase_name = gen_phase_name;
            gen_phase_name += "_efficiency";
            fill_unfold_hist(unfold_setup, UnfoldSpaceName::unfolded);
            gen_phase_name = original_gen_phase_name;

            if (_event.PassTrigger(p.triggers))
            {
                string original_gen_phase_name = gen_phase_name;
                gen_phase_name += "_trigger";
                fill_unfold_hist(unfold_setup, UnfoldSpaceName::unfolded);
                gen_phase_name = original_gen_phase_name;
            }
        }
    }
    else {
        return;
    }
}


void ISRUnfold::fill_unfold_reco_hist(ISRUnfoldSetUp* unfold_setup) {
    bool pass_reco = unfold_setup->pass_reco_cuts(reco_dipt, reco_dimass);
    if (pass_reco){
        fill_unfold_hist(unfold_setup, UnfoldSpaceName::folded);
    }
    else {
        return;
    }
}

void ISRUnfold::fill_unfold_hist(ISRUnfoldSetUp* unfold_setup, const UnfoldSpaceName mode){
    map<TString,double>* weights;
    if (mode == UnfoldSpaceName::unfolded) {
        weights = &gen_weights;
    } else {
        weights = &reco_weights;
    }
    for (const auto& [suffix,weight]:*weights){
        fill_unfold_hist(unfold_setup, suffix, weight, mode);
    }
}

 void ISRUnfold::fill_unfold_hist(ISRUnfoldSetUp* unfold_setup, TString suf, Double_t weight, const UnfoldSpaceName mode){
    // bin_name: [dipt-dimass]_[unfolded_fine_O-window_v1_UO]
    // "[" + dipt-dimass + "]_[" + level + "__unfolded_fine_O-window_v1_UO + "]" 
    if (unfold_setup->bin_turned_off(mode)){
        return;
    } 
    
    bool is_2d = unfold_setup->is_2d_unfold();
    string bin_name = unfold_setup->get_bin_name(mode);
    // check if this bin is turned off
    // then used return
    string phase_name = get_phase_name(mode);

    if (is_2d){
        // make hist name
        vector<string> matchs = get_reg_matchs(bin_name);
        string hname = string(p.prefix) + string(p.hprefix) + "[tunfold-hist]_" + "[" + matchs.at(0) + "]_[" + 
            phase_name + "__" + matchs.at(1) + "]" + string(p.suffix) + string(suf);

        // fill hist
        TUnfoldBinning* bin = map_unfold_2d_bins[bin_name];
        TH1D *this_hist = GetHist1D(hname);
        if (!this_hist){
            this_hist = (TH1D*) bin->CreateHistogram(hname.data());
            this_hist->SetDirectory(NULL);
            maphist_TH1D[hname] = this_hist;
        }

        int index = get_bin_index(bin, mode);
        this_hist->Fill(index, weight);
    }
    else {
        // dipt_[reco__fine_O]_[cuts]_dimass_55.0to64.0 
        string var_name = unfold_setup->get_var_name();
        string second_var_range;
        if (var_name == "dipt")
            second_var_range = unfold_setup->get_passed_winow_name(mode);
        else // dimass
            second_var_range = unfold_setup->get_dipt_range();
        string hname = string(p.prefix) + string(p.hprefix) + var_name + 
            "_[" + phase_name + "__" + unfold_setup->get_raw_bin_name(mode) + "]_" + 
            unfold_setup->get_second_axis_var_name() + "_" + second_var_range + string(p.suffix) + string(suf);
        // check here
        FillHist(hname, get_value(mode, var_name), weight, map_unfold_1d_bins[bin_name].size()-1, map_unfold_1d_bins[bin_name].data());
    }
}
void ISRUnfold::WriteHist(){

    SMPAnalyzerCore::WriteHist();
    outfile->cd();
    // loop over bin definition
    if (write_bins==true){
        for (std::map<TString,TUnfoldBinning*>::iterator mapit = map_unfold_2d_bins.begin(); mapit!=map_unfold_2d_bins.end(); mapit++){
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
