#include "ISRUnfold.h"

ISRUnfold::ISRUnfold(){

}

ISRUnfold::~ISRUnfold(){

}

void ISRUnfold::fill_unfold_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> weights, const TUnfoldParameter& par, const int mode){

    TLorentzVector dilepton = (*l0) + (*l1);
    double dimass = dilepton.M();
    double dipt = dilepton.Pt();

    int index{-1};
    TUnfoldBinning* bin_pointer = nullptr;

    std::map<TString, tuple<TUnfoldBinning*, TUnfoldBinning*>>::iterator mapit = map_tunfoldbins.find(par.bin_name);
    if(mapit != map_tunfoldbins.end()){  // bin definiton exists

        if(mode == 0) bin_pointer = get<0>(map_tunfoldbins[par.bin_name]);
        else if(mode == 1) bin_pointer = get<1>(map_tunfoldbins[par.bin_name]);
    }
    else{

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

        if(mode == 0) bin_pointer = get<0>(map_tunfoldbins[par.bin_name]);
        else if(mode == 1) bin_pointer = get<1>(map_tunfoldbins[par.bin_name]);
    }

    string var_name;  
    if(par.is_2D == true){
        if(bin_pointer->GetDistributionAxisLabel(0) == "dipt"){
            index = bin_pointer->GetGlobalBinNumber(dipt, dimass);
            var_name = "pt_mass";
        }
        else{
            index = bin_pointer->GetGlobalBinNumber(dimass, dipt);
            var_name = "mass_pt";
        }
    }

    string bin_prefix;
    if(mode == 0) bin_prefix = "_smeared";
    else if(mode == 1) bin_prefix = "_truth";
    // check dimension of bin   
    fill_unfold_hist(channelname+"/"+pre+par.bin_name + "_" + var_name+bin_prefix, index, weights, bin_pointer);
}

void ISRUnfold::fill_unfold_hist(TString histname, Double_t value, map<TString,double> weights, TUnfoldBinning* bin_pointer){

    for(const auto& [suffix,weight]:weights){
        fill_unfold_hist(histname+suffix,value,weight, bin_pointer);
    }
}

void ISRUnfold::fill_unfold_hist(TString hname, Double_t value, Double_t weight, TUnfoldBinning* bin_pointer){
    TH1D *this_hist = GetHist1D(hname);
    if( !this_hist ){
        this_hist = (TH1D*) bin_pointer->CreateHistogram(hname);

        this_hist->SetDirectory(NULL);
        maphist_TH1D[hname] = this_hist;
  }
  this_hist->Fill(value, weight);
}


