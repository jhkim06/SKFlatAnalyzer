#include "ISRUnfold.h"

ISRUnfold::ISRUnfold(){
}

ISRUnfold::~ISRUnfold(){
}

bool ISRUnfold::is_same_bin(TString channelname, TString pre,
                            Double_t variable1, Double_t variable2, const TUnfoldParameter& par, const TUnfold_Bin mode)
{

    TUnfoldBinning* temp_bin = new TUnfoldBinning("temp_bin");
    temp_bin->AddAxis(par.first_axis_var_name, par.n_first_axis_unfolded, par.first_axis_unfolded.data(), par.use_first_axis_uf, par.use_first_axis_of);
    
    TH1D* temp_hist = (TH1D*) temp_bin->CreateHistogram("temp_hist", true);

    int bin_index1 = temp_hist->FindBin(variable1);
    int bin_index2 = temp_hist->FindBin(variable2);


    if (bin_index1==bin_index2 && bin_index1 > 0 && bin_index1 <= temp_hist->GetNbinsX())
    {
        delete temp_bin;
        delete temp_hist;

        return true;
    }
    else{
        delete temp_bin;
        delete temp_hist;

        return false;
    }

}

void ISRUnfold::fill_unfold_hists(TString channelname, TString pre, TString suf,
                                  Particle* l0, Particle* l1, map<TString,double> weights, const TUnfoldParameter& par, const TUnfold_Bin mode){

    TLorentzVector dilepton = (*l0) + (*l1);
    double dimass = dilepton.M();
    double dipt = dilepton.Pt();

    double index{-1};
    TUnfoldBinning* bin_pointer = nullptr;

    // get desired bin definition
    string full_bin_name = (string)channelname+(string)pre+par.bin_name;
    std::map<TString, tuple<TUnfoldBinning*, TUnfoldBinning*>>::iterator mapit = map_tunfoldbins.find(full_bin_name);
    
    if (mapit != map_tunfoldbins.end()){  // bin definiton exists
        if (mode == TUnfold_Bin::folded_bin) bin_pointer = get<0>(map_tunfoldbins[full_bin_name]);
        else if (mode == TUnfold_Bin::unfolded_bin) bin_pointer = get<1>(map_tunfoldbins[full_bin_name]);
    }
    else{
        // create bin definition TODO make function to create bin definition
        TUnfoldBinning* unfolded_bin = new TUnfoldBinning("unfolded");
        TUnfoldBinning* folded_bin = new TUnfoldBinning("folded");

        unfolded_bin->AddAxis(par.first_axis_var_name, par.n_first_axis_unfolded, par.first_axis_unfolded.data(), par.use_first_axis_uf, par.use_first_axis_of);
        folded_bin  ->AddAxis(par.first_axis_var_name, par.n_first_axis_folded, par.first_axis_folded.data(), par.use_first_axis_uf, par.use_first_axis_of);
        
        // for 2D, add the second axis
        if(par.is_2D){
        unfolded_bin->AddAxis(par.second_axis_var_name, par.n_second_axis_unfolded, par.second_axis_unfolded.data(), par.use_second_axis_uf, par.use_second_axis_of);
        folded_bin->  AddAxis(par.second_axis_var_name, par.n_second_axis_folded, par.second_axis_folded.data(), par.use_second_axis_uf, par.use_second_axis_of);
        }
        
        map_tunfoldbins[full_bin_name] = make_tuple(folded_bin, unfolded_bin);

        if (mode == TUnfold_Bin::folded_bin) bin_pointer = get<0>(map_tunfoldbins[full_bin_name]);
        else if (mode == TUnfold_Bin::unfolded_bin) bin_pointer = get<1>(map_tunfoldbins[full_bin_name]);
    }
    
    if(par.is_2D == true){
        if (bin_pointer->GetDistributionAxisLabel(0) == "dipt"){
            index = bin_pointer->GetGlobalBinNumber(dipt, dimass);
        }
        else {
            index = bin_pointer->GetGlobalBinNumber(dimass, dipt);
        }
    }
    else{
        // 1D binning
        const string pt_str = "dipt";
        const string mass_str = "dimass";

        if (full_bin_name.find(pt_str)!=string::npos){
            index = dipt;
        }
        else if (full_bin_name.find(mass_str)!=string::npos){
            index = dimass;
        }
        else {
            cout <<"ISRUnfold::fill_unfold_hists check bin definition." << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    string bin_prefix;
    if (mode == TUnfold_Bin::folded_bin) bin_prefix = "folded";
    else if (mode == TUnfold_Bin::unfolded_bin) bin_prefix = "unfolded";
    // check dimension of bin   
    fill_unfold_hist(channelname+pre+par.bin_name+"_"+bin_prefix+suf, index, weights, bin_pointer, par.is_2D);
}

void ISRUnfold::fill_unfold_response_matrixs(TString channelname, TString pre, TString suf,
                                             Particle* l0, Particle* l1, Particle* unfolded_l0, Particle* unfolded_l1, map<TString,double> reco_weights, map<TString,double> gen_weights, const TUnfoldParameter& par)
{

    TLorentzVector dilepton_folded=(*l0)+(*l1);
    double dimass_folded=dilepton_folded.M();
    double dipt_folded=dilepton_folded.Pt();
    
    TLorentzVector dilepton_unfolded=(*unfolded_l0)+(*unfolded_l1);
    double dimass_unfolded = dilepton_unfolded.M();
    double dipt_unfolded = dilepton_unfolded.Pt();
    
    TUnfoldBinning* bin_pointer_folded = nullptr;
    TUnfoldBinning* bin_pointer_unfolded = nullptr;
    
    string full_bin_name = (string)channelname+(string)pre+par.bin_name;
    std::map<TString, tuple<TUnfoldBinning*, TUnfoldBinning*>>::iterator mapit = map_tunfoldbins.find(full_bin_name);
    
    if (mapit != map_tunfoldbins.end()){  // bin definiton exists
        bin_pointer_folded=get<0>(map_tunfoldbins[full_bin_name]);
        bin_pointer_unfolded=get<1>(map_tunfoldbins[full_bin_name]);
    }
    else {
        // create bin definition
        TUnfoldBinning* unfolded_bin = new TUnfoldBinning("unfolded"); // TODO get bin name from bin parameter
        TUnfoldBinning* folded_bin = new TUnfoldBinning("folded");

        unfolded_bin->AddAxis(par.first_axis_var_name, par.n_first_axis_unfolded, par.first_axis_unfolded.data(), par.use_first_axis_uf, par.use_first_axis_of);
        folded_bin->AddAxis(par.first_axis_var_name, par.n_first_axis_folded, par.first_axis_folded.data(), par.use_first_axis_uf, par.use_first_axis_of);
        
        if(par.is_2D){
        unfolded_bin->AddAxis(par.second_axis_var_name, par.n_second_axis_unfolded, par.second_axis_unfolded.data(), par.use_second_axis_uf, par.use_second_axis_of);
        folded_bin->AddAxis(par.second_axis_var_name, par.n_second_axis_folded, par.second_axis_folded.data(), par.use_second_axis_uf, par.use_second_axis_of);
        }

        map_tunfoldbins[full_bin_name] = make_tuple(folded_bin, unfolded_bin);
        
        bin_pointer_folded=get<0>(map_tunfoldbins[full_bin_name]);
        bin_pointer_unfolded=get<1>(map_tunfoldbins[full_bin_name]);
    }
    
    double index_folded{-1}, index_unfolded{-1};
    if (par.is_2D == true){
        if (bin_pointer_folded->GetDistributionAxisLabel(0)=="dipt"){
            index_folded=bin_pointer_folded->GetGlobalBinNumber(dipt_folded, dimass_folded);
            index_unfolded=bin_pointer_unfolded->GetGlobalBinNumber(dipt_unfolded, dimass_unfolded);
        }
        else {
            index_folded=bin_pointer_folded->GetGlobalBinNumber(dimass_folded, dipt_folded);
            index_unfolded=bin_pointer_unfolded->GetGlobalBinNumber(dimass_unfolded, dipt_unfolded);
        }
    }
    else {
        
        const string pt_str = "dipt";
        const string mass_str = "dimass";

        if (full_bin_name.find(pt_str)!=string::npos){

            index_folded = dipt_folded;
            index_unfolded = dipt_unfolded;
        }
        else if (full_bin_name.find(mass_str)!=string::npos){
            index_folded = dimass_folded;
            index_unfolded = dimass_unfolded;
        }
        else {
            cout <<"ISRUnfold::fill_unfold_hists check bin definition." << endl;
            exit(EXIT_FAILURE);
        }
    }

    string bin_prefix="responseM";
    fill_unfold_response_matrix(channelname+pre+par.bin_name+"_"+bin_prefix+suf, index_folded, index_unfolded,
                     reco_weights, gen_weights, bin_pointer_folded, bin_pointer_unfolded, par.is_2D);
}

void ISRUnfold::fill_unfold_hist(TString histname, Double_t value, map<TString,double> weights, TUnfoldBinning* bin_pointer, bool is_2D){

    for (const auto& [suffix,weight]:weights){
        fill_unfold_hist(histname+suffix, value, weight, bin_pointer, is_2D);
    }
}

void ISRUnfold::fill_unfold_response_matrix(TString hname, Double_t value_folded, Double_t value_unfolded, map<TString,double> reco_weights, map<TString,double> gen_weights, TUnfoldBinning* bin_pointer_folded, TUnfoldBinning* bin_pointer_unfolded, bool is_2D)
{
    for (const auto& [suffix,reco_weight]:reco_weights){
        fill_unfold_response_matrix(hname+suffix, value_folded, value_unfolded, reco_weight, gen_weights[suffix], bin_pointer_folded, bin_pointer_unfolded, is_2D);
    }
}

void ISRUnfold::fill_unfold_hist(TString hname, Double_t value, Double_t weight, TUnfoldBinning* bin_pointer, bool is_2D){
    TH1D *this_hist = GetHist1D(hname);
    if (!this_hist){

        this_hist = (TH1D*) bin_pointer->CreateHistogram(hname, !is_2D);

        this_hist->SetDirectory(NULL);
        maphist_TH1D[hname] = this_hist;
  }
  this_hist->Fill(value, weight);
}

void ISRUnfold::fill_unfold_response_matrix(TString hname, Double_t value_folded, Double_t value_unfolded, Double_t reco_weight, Double_t gen_weight, TUnfoldBinning* bin_pointer_folded, TUnfoldBinning* bin_pointer_unfolded, bool is_2D)
{
    TH2D *this_hist = GetHist2D(hname);

    if (!this_hist){
        this_hist = (TH2D*) TUnfoldBinning::CreateHistogramOfMigrations(bin_pointer_unfolded, bin_pointer_folded, hname, !is_2D, !is_2D);

        this_hist->SetDirectory(NULL);
        maphist_TH2D[hname] = this_hist;
  }
    
  this_hist->Fill(value_unfolded, value_folded, reco_weight);
  if (is_2D) this_hist->Fill(value_unfolded, 0., gen_weight-reco_weight); // bin zero for 2D
  else this_hist->Fill(value_unfolded, -1., gen_weight-reco_weight); // bin zero for 1D
}

void ISRUnfold::WriteHist(){

    SMPAnalyzerCore::WriteHist();
    outfile->cd();
    // loop over bin definition

    if (write_bins==true){
        for (std::map<TString,tuple<TUnfoldBinning*, TUnfoldBinning*>>::iterator mapit = map_tunfoldbins.begin(); mapit!=map_tunfoldbins.end(); mapit++){
            
            TString this_fullname=mapit->first;
            TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
            TString this_suffix=this_fullname(0,this_fullname.Last('/'));
            TDirectory *dir = outfile->GetDirectory(this_suffix);
            if(!dir){
                outfile->mkdir(this_suffix);
            }
            outfile->cd(this_suffix);
            get<0>(mapit->second)->Write(this_name+"_folded_bin");
            get<1>(mapit->second)->Write(this_name+"_unfolded_bin");
            outfile->cd();
        }
    }
}
