#include "ISRUnfold.h"

ISRUnfold::ISRUnfold(){
    
}

ISRUnfold::~ISRUnfold(){
    
}
/*
TH1D* ISRUnfold::get_hist1D(TString histname){

    TH1D *h = NULL;
    std::map<TString, TH1D*>::iterator mapit = maphistunfold_TH1D.find(histname);
    if(mapit != maphistunfold_TH1D.end()) return mapit->second;

    return h;
}
*/
void ISRUnfold::fill_unfold_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> weights, bool use_coarse_hist){

    TLorentzVector dilepton = (*l0) + (*l1);
    double dimass = dilepton.M();
    double dipt = dilepton.Pt();
    
    int index{-1};
    
    if(is_2D){
        TUnfoldBinning* bin_pointer = nullptr;
        
        // select bin definition
        string bin_prefix;
        if(!use_coarse_hist){
            bin_pointer = fine_bin;
            bin_prefix = "_fine";
        }
        else{
            bin_pointer = coarse_bin;
            bin_prefix = "_coarse";
        }
       
        // get bin index to fill
        string var_name;
        if(bin_pointer->GetDistributionAxisLabel(0) == "dipt"){
            index = bin_pointer->GetGlobalBinNumber(dipt, dimass);
            var_name = "pt_mass";
        }
        else{
            index = bin_pointer->GetGlobalBinNumber(dimass, dipt);
            var_name = "mass_pt";
        }
        
        fill_unfold_hist(channelname+"/"+pre+var_name+bin_prefix, index, weights);
    }
    
    
}

void ISRUnfold::fill_unfold_hist(TString histname, Double_t value, map<TString,double> weights){
    
    // append coarse, or fine in the histogram name
    for(const auto& [suffix,weight]:weights){
        fill_unfold_hist(histname+suffix,value,weight);
    }
}

void ISRUnfold::fill_unfold_hist(TString hname, Double_t value, Double_t weight){
    TH1D *this_hist = GetHist1D(hname);
    if( !this_hist ){
        if(hname.Contains("coarse")){
            this_hist = (TH1D*) coarse_bin->CreateHistogram(hname);
        }
        else{
            this_hist = (TH1D*) fine_bin->CreateHistogram(hname);
        }
      
        this_hist->SetDirectory(NULL);
        maphist_TH1D[hname] = this_hist;
  }
  this_hist->Fill(value, weight);
}

    
