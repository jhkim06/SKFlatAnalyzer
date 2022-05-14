#include "ISRUnfold.h"

ISRUnfold::ISRUnfold(){

}

ISRUnfold::~ISRUnfold(){

}

void ISRUnfold::fill_unfold_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> weights, const TUnfoldParameter& par, const TUnfold_Bin mode){

    TLorentzVector dilepton = (*l0) + (*l1);
    double dimass = dilepton.M();
    double dipt = dilepton.Pt();

    int index{-1};
    TUnfoldBinning* bin_pointer = nullptr;

    // get desired bin definition
    std::map<TString, tuple<TUnfoldBinning*, TUnfoldBinning*>>::iterator mapit = map_tunfoldbins.find(par.bin_name);
    if (mapit != map_tunfoldbins.end()){  // bin definiton exists

        if (mode == TUnfold_Bin::smeared_bin) bin_pointer = get<0>(map_tunfoldbins[par.bin_name]);
        else if (mode == TUnfold_Bin::truth_bin) bin_pointer = get<1>(map_tunfoldbins[par.bin_name]);
    }
    else{
        // create bin definition TODO make function to create bin definition
        TUnfoldBinning* truth_bin = new TUnfoldBinning("truth");
        TUnfoldBinning* smeared_bin = new TUnfoldBinning("smeared");

        truth_bin->AddAxis(par.first_axis_name, par.n_first_axis_truth, par.first_axis_truth, par.use_first_axis_uf, par.use_first_axis_of);
        smeared_bin->AddAxis(par.first_axis_name, par.n_first_axis_smeared, par.first_axis_smeared, par.use_first_axis_uf, par.use_first_axis_of);
        
        if(par.is_2D){
        truth_bin->AddAxis(par.second_axis_name, par.n_second_axis_truth, par.second_axis_truth, par.use_second_axis_uf, par.use_second_axis_of);
        smeared_bin->AddAxis(par.second_axis_name, par.n_second_axis_smeared, par.second_axis_smeared, par.use_second_axis_uf, par.use_second_axis_of);
        }

        map_tunfoldbins[par.bin_name] = make_tuple(smeared_bin, truth_bin);

        if (mode == TUnfold_Bin::smeared_bin) bin_pointer = get<0>(map_tunfoldbins[par.bin_name]);
        else if (mode == TUnfold_Bin::truth_bin) bin_pointer = get<1>(map_tunfoldbins[par.bin_name]);
    }

    string var_name;  
    if (par.is_2D == true){
        if (bin_pointer->GetDistributionAxisLabel(0) == "dipt"){
            index = bin_pointer->GetGlobalBinNumber(dipt, dimass);
            var_name = "pt_mass";
        }
        else {
            index = bin_pointer->GetGlobalBinNumber(dimass, dipt);
            var_name = "mass_pt";
        }
    }

    string bin_prefix;
    if (mode == TUnfold_Bin::smeared_bin) bin_prefix = "_smeared";
    else if (mode == TUnfold_Bin::truth_bin) bin_prefix = "_truth";
    // check dimension of bin   
    fill_unfold_hist(channelname+"/"+pre+par.bin_name + "_" + var_name+bin_prefix, index, weights, bin_pointer);
}

void ISRUnfold::fill_unfold_response_matrixs(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, Particle* truth_l0, Particle* truth_l1, map<TString,double> reco_weights, map<TString,double> gen_weights, const TUnfoldParameter& par)
{

    TLorentzVector dilepton_smeared=(*l0)+(*l1);
    double dimass_smeared=dilepton_smeared.M();
    double dipt_smeared=dilepton_smeared.Pt();
    
    TLorentzVector dilepton_truth=(*truth_l0)+(*truth_l1);
    double dimass_truth = dilepton_truth.M();
    double dipt_truth = dilepton_truth.Pt();
    
    TUnfoldBinning* bin_pointer_smeared = nullptr;
    TUnfoldBinning* bin_pointer_truth = nullptr;
    
    std::map<TString, tuple<TUnfoldBinning*, TUnfoldBinning*>>::iterator mapit = map_tunfoldbins.find(par.bin_name);
    if (mapit != map_tunfoldbins.end()){  // bin definiton exists
        bin_pointer_smeared=get<0>(map_tunfoldbins[par.bin_name]);
        bin_pointer_truth=get<1>(map_tunfoldbins[par.bin_name]);
    }
    else {
        // create bin definition
        TUnfoldBinning* truth_bin = new TUnfoldBinning("truth");
        TUnfoldBinning* smeared_bin = new TUnfoldBinning("smeared");

        truth_bin->AddAxis(par.first_axis_name, par.n_first_axis_truth, par.first_axis_truth, par.use_first_axis_uf, par.use_first_axis_of);
        smeared_bin->AddAxis(par.first_axis_name, par.n_first_axis_smeared, par.first_axis_smeared, par.use_first_axis_uf, par.use_first_axis_of);
        
        if(par.is_2D){
        truth_bin->AddAxis(par.second_axis_name, par.n_second_axis_truth, par.second_axis_truth, par.use_second_axis_uf, par.use_second_axis_of);
        smeared_bin->AddAxis(par.second_axis_name, par.n_second_axis_smeared, par.second_axis_smeared, par.use_second_axis_uf, par.use_second_axis_of);
        }

        map_tunfoldbins[par.bin_name] = make_tuple(smeared_bin, truth_bin);
        
        bin_pointer_smeared=get<0>(map_tunfoldbins[par.bin_name]);
        bin_pointer_truth=get<1>(map_tunfoldbins[par.bin_name]);
    }
    
    int index_smeared{-1}, index_truth{-1};
    string var_name;
    if (par.is_2D == true){
        if (bin_pointer_smeared->GetDistributionAxisLabel(0)=="dipt"){
            index_smeared=bin_pointer_smeared->GetGlobalBinNumber(dipt_smeared, dimass_smeared);
            index_truth=bin_pointer_truth->GetGlobalBinNumber(dipt_truth, dimass_truth);
            var_name="pt_mass";
        }
        else {
            index_smeared=bin_pointer_smeared->GetGlobalBinNumber(dimass_smeared, dipt_smeared);
            index_truth=bin_pointer_truth->GetGlobalBinNumber(dimass_truth, dipt_truth);
            var_name = "mass_pt";
        }
    }
    string bin_prefix="res_matrix";
    fill_unfold_response_matrix(channelname+"/"+pre+par.bin_name+"_"+var_name+bin_prefix, index_smeared, index_truth,
                     reco_weights, gen_weights, bin_pointer_smeared, bin_pointer_truth);
}

void ISRUnfold::fill_unfold_hist(TString histname, Double_t value, map<TString,double> weights, TUnfoldBinning* bin_pointer){

    for (const auto& [suffix,weight]:weights){
        fill_unfold_hist(histname+suffix,value,weight, bin_pointer);
    }
}

void ISRUnfold::fill_unfold_response_matrix(TString hname, Int_t value_smeared, Int_t value_truth, map<TString,double> reco_weights, map<TString,double> gen_weights, TUnfoldBinning* bin_pointer_smeared, TUnfoldBinning* bin_pointer_truth)
{
    for (const auto& [suffix,reco_weight]:reco_weights){
        fill_unfold_response_matrix(hname+suffix, value_smeared, value_truth, reco_weight, gen_weights[suffix], bin_pointer_smeared, bin_pointer_truth);
    }
}

void ISRUnfold::fill_unfold_hist(TString hname, Double_t value, Double_t weight, TUnfoldBinning* bin_pointer){
    TH1D *this_hist = GetHist1D(hname);
    if (!this_hist){
        this_hist = (TH1D*) bin_pointer->CreateHistogram(hname);

        this_hist->SetDirectory(NULL);
        maphist_TH1D[hname] = this_hist;
  }
  this_hist->Fill(value, weight);
}

void ISRUnfold::fill_unfold_response_matrix(TString hname, Int_t value_smeared, Int_t value_truth, Double_t reco_weight, Double_t gen_weight, TUnfoldBinning* bin_pointer_smeared, TUnfoldBinning* bin_pointer_truth)
{
    TH2D *this_hist = GetHist2D(hname);

    if (!this_hist){
        this_hist = (TH2D*) TUnfoldBinning::CreateHistogramOfMigrations(bin_pointer_truth, bin_pointer_smeared, hname);

        this_hist->SetDirectory(NULL);
        maphist_TH2D[hname] = this_hist;
  }
    
  this_hist->Fill(value_truth, value_smeared, reco_weight*gen_weight);
  this_hist->Fill(value_truth, 0., (1-reco_weight)*gen_weight); // bin zero
    
}


