#include "ISRAnalyzer.h"

void ISRAnalyzer::initializeAnalyzer(){

    SMPAnalyzerCore::initializeAnalyzer(); // setup zpt roc z0
    ISRUnfold::initializeISRUnfold(job_number); // to write TUnfoldBinning to the first output root file

    // FIXME already in SMPAnalyzerCore::initializeAnalyzer()?
    vector<JetTagging::Parameters> jtps={
        JetTagging::Parameters(JetTagging::DeepJet,
                               JetTagging::Tight,
                               JetTagging::incl,
                               JetTagging::comb),

        JetTagging::Parameters(JetTagging::DeepJet,
                               JetTagging::Medium,
                               JetTagging::incl,
                               JetTagging::comb),

        JetTagging::Parameters(JetTagging::DeepJet,
                               JetTagging::Loose,
                               JetTagging::incl,
                               JetTagging::comb),
    };
    mcCorr->SetJetTaggingParameters(jtps);

    if (fChain->GetListOfFiles()->GetEntries()) {
        TString filename=fChain->GetListOfFiles()->At(0)->GetTitle();
        if(filename.Contains("SkimTree_")) IsSkimmed=true;
        else IsSkimmed=false;
    } else {
        cout<<"[ISRAnalyzer::initializeAnalyzer] no input file"<<endl;
        exit(EXIT_FAILURE);
    }

    IsNominalRun=!HasFlag("SYS")&&!HasFlag("PDFSYS");
    // TODO Flag for unfolding validation
}

void ISRAnalyzer::executeEvent(){
    // FIXME some events of DYJets has nan PDF weights. I don't know why...
    if(MCSample=="DYJets"&&!isnormal(weight_Scale->at(0))) return;

    // TODO effect of trigger matching?
    if(!IsDATA||DataStream.Contains("DoubleMuon")){
        //executeEventWithParameter(MakeParameter("mm", "TightID_TightIso"));
        executeEventWithParameter(MakeParameter("mm", "TightID_TightIso_b_veto"));
        if(IsNominalRun){
            executeEventWithParameter(MakeParameter("mm", "TightID_TightIso_b_veto",  "roccor_s2m0")); //Zpt
            executeEventWithParameter(MakeParameter("mm", "TightID_TightIso_b_veto",  "roccor_s3m0")); //Ewk
            executeEventWithParameter(MakeParameter("mm", "TightID_TightIso_b_veto",  "roccor_s4m0")); //deltaM
            executeEventWithParameter(MakeParameter("mm", "TightID_TightIso_b_veto",  "roccor_s5m0")); //Ewk2
			// stat
			for (int N = 0; N < 100; ++N) {
				TString rocParam = TString::Format("roccor_s1m%d", N);
				executeEventWithParameter(MakeParameter("mm","TightID_TightIso_b_veto",rocParam));
			}
        }
    }
    if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
        executeEventWithParameter(MakeParameter("ee", "TightID_b_veto"));
		if(IsNominalRun){
			executeEventWithParameter(MakeParameter("ee", "TightID_b_veto", "egammacor_s1m-1"));  // FIXME egamma_s1m-1 is not shared option actually
			executeEventWithParameter(MakeParameter("ee", "TightID_b_veto", "egammacor_s1m1"));
			executeEventWithParameter(MakeParameter("ee", "TightID_b_veto", "egammacor_s2m-1"));
			executeEventWithParameter(MakeParameter("ee", "TightID_b_veto", "egammacor_s2m1"));
		}
        // executeEventWithParameter(MakeParameter("ee", "MediumID_b_veto"));
        // executeEventWithParameter(MakeParameter("ee", "TightID_b_veto_remove_EE"));
        if (DataYear == 2018) {
            //executeEventWithParameter(MakeParameter("ee", "TightID_no_vetoHEM"));
            //executeEventWithParameter(MakeParameter("ee", "TightID_b_veto_no_vetoHEM"));
        }
    }
}

SMPAnalyzerCore::Parameter ISRAnalyzer::MakeParameter(TString key, TString option, TString option_shared){ //

    // CAUTION: avoid option name defined in SMPAnalyzer
    Parameter p = SMPAnalyzerCore::MakeParameter(key, option_shared);
    // Update SMPAnalyzerCore::MakeParameter
    // ELECTRON channel
    if (key == "ee"){
        TString electron_id_name;
        bool vetoHEM = true; // only for 2018
        double eta_cut = 2.5;

        // Lepton ID
        if (option.Contains("TightID"))
        {
            electron_id_name = "passTightID";
            p.SetElectronKeys("Electron_TightID",{"Ele23Leg1_TightID","Ele12Leg2_TightID"});
            if (GetEraShort()=="2016a"){
                p.k.DZSF="DZ_MediumID";  // FIXME
            }
        }
        else
        {
            electron_id_name = "passMediumID";
            if (!option.Contains("MediumID")) {
                option = "MediumID_" + option;
            }
        }
        if (option.Contains("no_vetoHEM")) {
            vetoHEM = false;
        }

        int set = 0;
        int member = 0;
		// build your regex once
		static TPRegexp re("egammacor_s([0-9]+)m(-1|1)");
		if ( re.Match(p.option) ) {
			// re.MatchS returns a new TObjArray* that you own
			TObjArray* arr = re.MatchS(p.option);
			if (arr && arr->GetEntries() >= 3) {
				// At(0) is the full match, At(1) is ([0-9]+), At(2) is (-1|1)
				TString sset   = ((TObjString*)arr->At(1))->GetString();
				TString smem   = ((TObjString*)arr->At(2))->GetString();
				set  = sset.Atoi();
				member  = smem.Atoi();
				p.suffix+="_"+((TObjString*)arr->At(0))->GetString(); 
			} else {
				// no groups? you might log or handle the error
			}
			delete arr;
		}
        p.SetLeptonEtaCut(eta_cut); // to be used later in the unfold histogram
        std::vector<Electron> electrons_ = ISRGetElectrons(electron_id_name, 0.0, eta_cut, vetoHEM, set , member);
        //std::vector<Electron> electrons_ = SMPAnalyzerCore::SMPGetElectrons("passTightID_SelQ", 0.0, eta_cut);

        if (option.Contains("remove_EE"))
        {
            std::vector<Electron> temp_electrons;
            for (unsigned int i = 0; i < electrons_.size(); i++) {
                if (RemoveElectron(electrons_.at(i)))
                    continue;
                temp_electrons.push_back(electrons_.at(i));
            }
            electrons_ = temp_electrons;
        }

        std::sort(electrons_.begin(), electrons_.end(), PtComparing);
        p.SetElectrons(electrons_);
    }
    // MUON channel
	/*
    if (key == "mm"){
        TString muon_id_iso_name;
        double eta_cut = 2.4;
        if (option.Contains("TightID_TightIso"))
        {
            muon_id_iso_name = "POGTightWithTightIso";
        }
        else {
            // default
            muon_id_iso_name = "POGTightWithTightIso";
            if (!option.Contains("TightID_TightIso")) {
                option = "TightID_TightIso_" + option;
            }
        }
        int sys = 0;
        int set = 0;
        int member = 0;
        std::vector<Muon> muons_ = MuonMomentumCorrection(SMPGetMuons(muon_id_iso_name, 0.0, eta_cut), sys, set, member);
        // change muon key
        p.SetMuons(muons_);
        p.SetMuonKeys("Muon_TightID_TightIso","",{"Mu17Leg1_TightID_TightIso","Mu8Leg2_TightID_TightIso"});
        // TODO NEED to CHANGE Trigger SF
        p.SetLeptonEtaCut(eta_cut);
    }
	*/
	if (key == "mm") {
		// 1) ensure option starts with TightID_TightIso_
		const TString prefix = "TightID_TightIso_";
		if (!option.BeginsWith(prefix)) {
			option = prefix + option;
		}

		// 2) the Muon ID+Iso working point is always the same
		const TString muon_id_iso_name = "POGTightWithTightIso";

		// 3) your usual parameters
		const double eta_cut = 2.4;
		const int sys    = 0;
		int set    = 0;
		int member = 0;

		// for systematic of rochester correction
		// p.suffix is  set at SMPAnalyzerCore::MakeParameter
  	    if(p.option.Contains(TRegexp("roccor_s[0-9]+m[0-9]+"))){
  	      	TObjArray* array=TPRegexp("roccor_s([0-9]+)m([0-9]+)").MatchS(p.option);
  	      	TString sset=((TObjString*)array->At(1))->GetString();
  	      	TString smem=((TObjString*)array->At(2))->GetString();
  	      	set=sset.Atoi();
  	      	member=smem.Atoi();
  	      	delete array;
  	    }

		// 4) grab & correct the muons
		auto raw_muons = SMPGetMuons(muon_id_iso_name, /*ptMin=*/0.0, eta_cut);
		auto muons     = MuonMomentumCorrection(raw_muons, sys, set, member);

		// 5) hand them off to your plotter 'p'
		p.SetMuons(muons);

		// 6) define your muon-trigger keys only once
		std::vector<TString> trigKeys = {
			"Mu17Leg1_TightID_TightIso",
			"Mu8Leg2_TightID_TightIso"
		};
		p.SetMuonKeys("Muon_TightID_TightIso", /*suffix=*/"", trigKeys);

		// 7) final bit of book-keeping
		p.SetLeptonEtaCut(eta_cut);
	}
    p.prefix += option + "/";

    p.weightbit = 0;
    if(IsNominalRun) p.weightbit|=NominalWeight;
    if(HasFlag("SYS")&&!IsDATA&&(p.channel=="ee"||p.channel=="mm"))
        p.weightbit|=SystematicWeight|EfficiencyWeight;
    if(HasFlag("PDFSYS")&&!IsDATA&&(p.channel=="ee"||p.channel=="mm"))
        p.weightbit|=PDFWeight;

    //if(HasFlag("nbjet")) p.prefix+="nbjet/";
    //else if(HasFlag("0bjet")) p.prefix+="0bjet/";
    //if(HasFlag("highmet")) p.prefix+="highmet/";

    p.option += option;  // "b_veto" used in PassSelection
    return p;
}

std::vector<Electron> ISRAnalyzer::ISRGetElectrons(TString id, double ptmin, double fetamax, bool vetoHEM, int set, int member) const {
  std::vector<Electron> electrons = GetAllElectrons();
  std::vector<Electron> out;
  for(unsigned int i=0; i<electrons.size(); i++){
	// set==1 scale, set==2 resolution
    // member==1 up, member==-1 down
	if(set==1){
		double scale = electrons.at(i).EnShift(member);
		electrons.at(i).SetPtEtaPhiE( scale*electrons.at(i).Pt(), electrons.at(i).Eta(), electrons.at(i).Phi(), scale*electrons.at(i).E());
	}
	else if(set==2){
		double res = electrons.at(i).ResShift(member);
		electrons.at(i).SetPtEtaPhiE( res*electrons.at(i).Pt(), electrons.at(i).Eta(), electrons.at(i).Phi(), res*electrons.at(i).E());
	}

    if(!( electrons.at(i).Pt()>ptmin )){
      //cout << "Fail Pt : pt = " << electrons.at(i).Pt() << ", cut = " << ptmin << endl;
      continue;
    }
    if(!( fabs(electrons.at(i).scEta())<fetamax )){
      //cout << "Fail Eta : eta = " << fabs(electrons.at(i).scEta()) << ", cut = " << fetamax << endl;
      continue;
    }
    if(!( electrons.at(i).PassID(id) )){
      //cout << "Fail ID" << endl;
      continue;
    }
    if(vetoHEM){
      if ( FindHEMElectron (electrons.at(i)) ){
        continue;
      }
    }
    out.push_back( electrons.at(i) );
  }
  return out;
}

bool ISRAnalyzer::RemoveElectron(Electron electron){
    if (fabs(electron.Eta()) > 2.0){
         return true;
    }
    return false;
}

bool ISRAnalyzer::RemoveMuon(Muon muon, double phi_min, double phi_max,
                             bool for_endcap, int eta_sign){
    // endcap_sign == 0
    // endcap_sign == 1
    // endcap_sign == -1
    bool pass_sign = true;
    if (eta_sign == 1) {
        if (muon.Eta() < 0) pass_sign = false;
    }
    if (eta_sign == -1) {
        if (muon.Eta() > 0) pass_sign = false;
    }

    bool is_endcap = false;
    if (fabs(muon.Eta()) > 0.9) is_endcap = true;

    if ( pass_sign && for_endcap == is_endcap && (phi_min <= muon.Phi() && muon.Phi() < phi_max))
        return true;
    else
        return false;
}

// called in SMPAnalyzerCore::executeEventWithParameter()
bool ISRAnalyzer::PassSelection(Parameter& p){
    if(!SMPAnalyzerCore::PassSelection(p)) return false;

    int n_bjet=0;
    std::vector<Jet> jets=GetJets("tightLepVeto", 30, 2.4);
    std::sort(jets.begin(), jets.end(), PtComparing);

    JetTagging::Parameters jtp;
    JetTagging::WP wp = JetTagging::Tight;
    jtp = JetTagging::Parameters(JetTagging::DeepJet,
                                 wp,
                                 JetTagging::incl,
                                 JetTagging::comb);

    for(const auto& jet:jets){
        if (jet.GetTaggerResult(jtp.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtp.j_Tagger,
                                                                              jtp.j_WP)) {
            n_bjet++;
        }
    }

    p.intmap["nbjet"]=n_bjet;
    p.doublemap["btagSF"]=mcCorr->GetBTaggingReweight_1a(jets, jtp);
    p.doublemap["btagSF_hup"]=mcCorr->GetBTaggingReweight_1a(jets, jtp, "SystUpHTag");
    p.doublemap["btagSF_hdown"]=mcCorr->GetBTaggingReweight_1a(jets, jtp, "SystDownHTag");
    p.doublemap["btagSF_lup"]=mcCorr->GetBTaggingReweight_1a(jets, jtp," SystUpLTag");
    p.doublemap["btagSF_ldown"]=mcCorr->GetBTaggingReweight_1a(jets, jtp, "SystDownLTag");

    if(n_bjet > 0 && p.option.Contains("b_veto")) return false; // b-jet veto
    //if(p.prefix.Contains("nbjet")&&!n_bjet) return false;
    //if(p.prefix.Contains("0bjet")&&n_bjet) return false;

    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCut", p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight);
    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCutSF", p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.doublemap["btagSF"]);
    return true;
}

void ISRAnalyzer::executeEventGen(Parameter& p){

    if(IsDYSample){

        DY_index = -1;
        gen_p0_isr = Gen();
        gen_p1_isr = Gen();
        gen_l0_bare = Gen();
        gen_l1_bare = Gen();
        gen_leptons.clear();
        gen_photons.clear();

        if(abs(lhe_l0.ID())!=15&&abs(lhe_l1.ID())!=15){
            if( (abs(lhe_l0.ID())==11&&abs(lhe_l1.ID())==11) || (!lhes.size()&&abs(gen_l0.PID())==11&&abs(gen_l1.PID())==11) ){
                if (p.channel != "ee") return;
            } else if( (abs(lhe_l0.ID())==13&&abs(lhe_l1.ID())==13) || (!lhes.size()&&abs(gen_l0.PID())==13&&abs(gen_l1.PID())==13) ){
                if (p.channel != "mm") return;
            } else {
                exit(EXIT_FAILURE);
            }

            get_DY_gen_particles();
            Gen gen_l0_isr;
            Gen gen_l1_isr;
			get_DY_dressed_lepton_pair(gen_l0_isr, gen_l1_isr, DressedMode::MotherMatch);

            //DY_index = get_DY_gen_particles(gens, gen_isr_parton0, gen_isr_parton1,
            //                                bare_l0, bare_l1);
            //DY_index = get_DY_gen_particles(gens, gen_isr_parton0, gen_isr_parton1,
            //                                dressedp1_l0, dressedp1_l1);
            // SMPAnalyzerCore::GetAFBGenParticles(gens, gen_isr_parton0, gen_isr_parton1, gen_isr_l0, gen_isr_l1, 3);
            // cout << "pre FSR mass: " << (gen_isr_l0 + gen_isr_l1).M() << " post FSR mass: " << (bare_l0 + bare_l1).M() << endl;

            Parameter pgen=p;
            ResetRecoWeights(pgen);
            EvalWeights(pgen);
            set_gen_weights(pgen.weightmap);

            set_base_parameter(pgen, false);
            set_gen_leptons(gen_l0_isr, gen_l1_isr);
            set_phase_name(UnfoldSpaceName::unfolded, "gen_dressed");
            fill_unfold_gen_hists();

            double preFSR_dimass = (gen_l0_isr + gen_l1_isr).M();
            double preFSR_dipt = (gen_l0_isr + gen_l1_isr).Pt();

            // loop over mass windows
			double dipt_cut = 100.0;
			const auto& vec = bins.at("dimass").at("coarse");
			// loop over it
			for (size_t i = 0;  i + 1 < vec.size();  ++i) {
				double low  = vec[i];
				double high = vec[i+1];
				if(preFSR_dimass>=low && preFSR_dimass<=high){
					TString hname = TString::Format("dipt_gen_dressed_acceptance_dipt_%.1fto%.1f_dimass_%.1fto%.1f", 0.0, dipt_cut, low, high);
					FillHist(p.prefix + hname, preFSR_dipt, pgen.weightmap, int(dipt_cut-0), 0, dipt_cut);
				}
			}

			const auto& mass_window_vec = windows.at("dimass").at("window_v1");
			for (size_t i = 0;  i + 1 < mass_window_vec.size();  ++i) {
				double low  = mass_window_vec[i];
				double high = mass_window_vec[i+1];
				if(preFSR_dimass>=low && preFSR_dimass<=high && preFSR_dipt<dipt_cut){
					TString hname = TString::Format("dimass_gen_dressed_acceptance_dipt_%.1fto%.1f_dimass_%.1fto%.1f", 0.0, dipt_cut, low, high);
					FillHist(p.prefix + hname, preFSR_dimass, pgen.weightmap, int(high-low), low, high);
					hname = TString::Format("dipt_gen_dressed_acceptance_dipt_%.1fto%.1f_dimass_%.1fto%.1f", 0.0, dipt_cut, low, high);
					FillHist(p.prefix + hname, preFSR_dipt, pgen.weightmap, int(dipt_cut-0), 0, dipt_cut);
				}
			}

            //set_gen_leptons(bare_l0, bare_l1);
            //set_phase_name(UnfoldSpaceName::unfolded, "gen_bare");
            //fill_unfold_gen_hists();

            //set_gen_leptons(dressedp1_l0, dressedp1_l1);
            //set_phase_name(UnfoldSpaceName::unfolded, "gen_dressedp1");
            //fill_unfold_gen_hists();

            // response matrix from bare to dressed
            // set reco pt as bare lepton
            // set gen pt as dressed lepton
            // set reco and gen weight as genweight
            //set_reco_leptons(bare_l0, bare_l1);
            //set_gen_leptons(gen_isr_l0, gen_isr_l1);
            //set_reco_weights(pgen.weightmap);
            //set_phase_name(UnfoldSpaceName::folded, "gen_bare");
            //set_phase_name(UnfoldSpaceName::unfolded, "gen_dressed");
            //fill_unfold_matrixs(false);  // remove kinematic cut reqirement
            //fill_unfold_fake_hists(false);  // remove kinematic cut requirement

            //set_reco_leptons(dressedp1_l0, dressedp1_l1);
            //set_gen_leptons(gen_isr_l0, gen_isr_l1);
            //set_reco_weights(pgen.weightmap);
            //set_phase_name(UnfoldSpaceName::folded, "gen_dressedp1");
            //set_phase_name(UnfoldSpaceName::unfolded, "gen_dressed");
            //fill_unfold_matrixs(false);  // remove kinematic cut reqirement
            //fill_unfold_fake_hists(false);  // remove kinematic cut requirement

            /*
            // for additional FSR info
            double LHE_dimass = (lhe_l0 + lhe_l1).M();
            double postFSR_dimass = (bare_l0 + bare_l1).M();

            // dimass
            FillHist(p.prefix + "LHE_dimass", LHE_dimass, pgen.weightmap[""], bins.at("dimass").at("coarse").size()-1, bins.at("dimass").at("coarse").data());  // parton
            FillHist(p.prefix + "preFSR_dimass", preFSR_dimass, pgen.weightmap[""], bins.at("dimass").at("coarse").size()-1, bins.at("dimass").at("coarse").data());  // dressed
            FillHist(p.prefix + "postFSR_dimass", postFSR_dimass, pgen.weightmap[""], bins.at("dimass").at("coarse").size()-1, bins.at("dimass").at("coarse").data());  // bare

            // 2D mass vs delta R
            vector<double> deltaR_bins;
            for (int i=0; i <= 500; i++){
                deltaR_bins.push_back(i*0.02);
            }
            // loop over photons
            for (const auto photon: added_photons) {
                double deltaR0 = bare_l0.DeltaR(*photon);
                double deltaR1 = bare_l1.DeltaR(*photon);
                double deltaR =  deltaR0<deltaR1 ? deltaR0 : deltaR1;
                FillHist(p.prefix+"deltaR_preFSR_dimass", preFSR_dimass, deltaR, pgen.weightmap[""],
                        bins.at("dimass").at("coarse").size()-1, bins.at("dimass").at("coarse").data(),
                        deltaR_bins.size()-1, deltaR_bins.data());
                // photon momentum
                FillHist(p.prefix + "FSR_photon_momentum", photon->Pt(), pgen.weightmap[""], 500, 0., 500.);
                FillHist(p.prefix + "FSR_photon_deltaR", deltaR, pgen.weightmap[""], deltaR_bins.size()-1, deltaR_bins.data());
            }
            // pt mass matrix
            FillHist(p.prefix + "preFSR_mass_pt_matrix", preFSR_dimass, preFSR_dipt, pgen.weightmap[""], 2000, 0., 2000., 1000, 0., 1000.);
            */
        }
    }
}

void ISRAnalyzer::executeEventWithParameter(Parameter& p){
    // draw cuflow according to the lepton selection defined in p
    p.intmap["nbjet"]=1.0;
    p.doublemap["btagSF"]=1.0;
    p.doublemap["btagSF_hup"]=1.0;
    p.doublemap["btagSF_hdown"]=1.0;
    p.doublemap["btagSF_lup"]=1.0;
    p.doublemap["btagSF_ldown"]=1.0;


    executeEventGen(p);
    SMPAnalyzerCore::executeEventWithParameter(p);
}

void ISRAnalyzer::EvalWeights(Parameter& p){

    if(p.weightbit&NominalWeight){
        p.weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
    }
    if(p.weightbit&SystematicWeight){
        if(!IsDATA){
            p.weightmap["_zptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
            p.weightmap["_noPUweight"]=p.w.lumiweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
            p.weightmap["_PUweight_up"]=p.w.lumiweight*p.w.PUweight_up*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
            p.weightmap["_PUweight_down"]=p.w.lumiweight*p.w.PUweight_down*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;

            p.weightmap["_noprefireweight"]=p.w.lumiweight*p.w.PUweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
            p.weightmap["_prefireweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_up*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
            p.weightmap["_prefireweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_down*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;

            p.weightmap["_noz0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
            p.weightmap["_noweakweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;

            p.weightmap["_nobtagSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.w.topptweight;
            p.weightmap["_btagSF_hup"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF_hup"]*p.w.topptweight;
            p.weightmap["_btagSF_hdown"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF_hdown"]*p.w.topptweight;
            p.weightmap["_btagSF_lup"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF_lup"]*p.w.topptweight;
            p.weightmap["_btagSF_ldown"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF_ldown"]*p.w.topptweight;

            for(int j=0,nj=fEff->nreplica;j<nj;j++){
                double electronRECOSF=p.w.electronRECOSF_sys.size() ? p.w.electronRECOSF_sys[0][j] : 1.;
                double electronIDSF=p.w.electronIDSF_sys.size() ? p.w.electronIDSF_sys[0][j] : 1.;
                double muonIDSF=p.w.muonIDSF_sys.size() ? p.w.muonIDSF_sys[0][j] : 1.;
                double triggerSF=p.w.triggerSF_sys.size() ? p.w.triggerSF_sys[0][j] : 1.;
                p.weightmap[Form("_efficiencySF_stat%d",j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*electronRECOSF*electronIDSF*muonIDSF*p.w.muonISOSF*triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
            }
            for(int i=1,ni=p.w.electronRECOSF_sys.size();i<ni;i++){
                for(int j=0,nj=p.w.electronRECOSF_sys[i].size();j<nj;j++){
                    p.weightmap[Form("_electronRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF_sys[i][j]*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
                }
            }
            for(int i=1,ni=p.w.electronIDSF_sys.size();i<ni;i++){
                for(int j=0,nj=p.w.electronIDSF_sys[i].size();j<nj;j++){
                    p.weightmap[Form("_electronIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF_sys[i][j]*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
                }
            }
            for(int i=1,ni=p.w.muonIDSF_sys.size();i<ni;i++){
                for(int j=0,nj=p.w.muonIDSF_sys[i].size();j<nj;j++){
                    p.weightmap[Form("_muonIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF_sys[i][j]*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*p.w.topptweight;
                }
            }
            for(int i=1,ni=p.w.triggerSF_sys.size();i<ni;i++){
                for(int j=0,nj=p.w.triggerSF_sys[i].size();j<nj;j++){
                    p.weightmap[Form("_triggerSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_sys[i][j]*p.doublemap["btagSF"]*p.w.topptweight;
                }
            }
        }
    }
    if(p.weightbit&PDFWeight){
        for(unsigned int i=0;i<weight_Scale->size();i++){
            p.weightmap[Form("_scalevariation%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_Scale->at(i)*p.w.topptweight;
        }
        for(unsigned int i=0;i<weight_PDF->size();i++){
            p.weightmap[Form("_pdf%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_PDF->at(i)*p.w.topptweight;
        }
        if(weight_AlphaS->size()==2){
            p.weightmap["_alphaS_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_AlphaS->at(0)*p.w.topptweight;
            p.weightmap["_alphaS_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_AlphaS->at(1)*p.w.topptweight;
        }

        if(MCSample.Contains("MiNNLO")){
            p.weightmap["_sthw2_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_sthw2->at(0)*p.w.topptweight;
            p.weightmap["_sthw2_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_sthw2->at(2)*p.w.topptweight;
            p.weightmap["_largeptscales"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_largeptscales->at(0)*p.w.topptweight;
            p.weightmap["_q0_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_q0->at(0)*p.w.topptweight;
            p.weightmap["_q0_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_q0->at(2)*p.w.topptweight;
        }
    }
    return;
}

void ISRAnalyzer::ResetRecoWeights(Parameter& p){
    p.w.prefireweight=1.; p.w.prefireweight_up=1.; p.w.prefireweight_down=1.;
    p.w.z0weight=1.;
    p.w.electronRECOSF=1.;
    p.w.electronRECOSF_sys=fEff->GetStructure(p.k.electronRECOSF);
    p.w.electronIDSF=1.;
    p.w.electronIDSF_sys=fEff->GetStructure(p.k.electronIDSF);
    p.w.muonIDSF=1.;
    p.w.muonIDSF_sys=fEff->GetStructure(p.k.muonIDSF);
    p.w.muonISOSF=1.;
    p.w.muonISOSF_sys=fEff->GetStructure(p.k.muonISOSF);
    p.w.triggerSF=1.;
    p.w.triggerSF_sys=fEff->GetStructure(p.k.triggerSF[0]);
    p.doublemap["btagSF"]=1.;
}

void ISRAnalyzer::FillHists(Parameter& p){

    set_base_parameter(p);
    set_phase_name(UnfoldSpaceName::folded, "reco");

    if (IsDYSample && p.hprefix!="tau_"){
        //DY_index = get_DY_gen_particles(gens,
        //        gen_isr_parton0, gen_isr_parton1,
        //        dressed_p1_l0, dressed_p1_l1,
        //        DressedDRp1, added_photons_p1);

        //DY_index = get_DY_gen_particles(gens,
        //        gen_isr_parton0, gen_isr_parton1,
        //        bare_l0, bare_l1);

        Gen gen_l0_isr;
        Gen gen_l1_isr;
		get_DY_dressed_lepton_pair(gen_l0_isr, gen_l1_isr, DressedMode::MotherMatch);
        set_gen_leptons(gen_l0_isr, gen_l1_isr);
        set_phase_name(UnfoldSpaceName::unfolded, "gen_dressed");
        Parameter pgen=p;
        ResetRecoWeights(pgen);
        EvalWeights(pgen);
        set_gen_weights(pgen.weightmap);

        fill_unfold_matrixs();
        fill_unfold_fake_hists();

        //set_gen_leptons(bare_l0, bare_l1);
        //set_phase_name(UnfoldSpaceName::unfolded, "gen_bare");

        //fill_unfold_matrixs();
        //fill_unfold_fake_hists();

        //set_gen_leptons(dressed_p1_l0, dressed_p1_l1);
        //set_phase_name(UnfoldSpaceName::unfolded, "gen_dressedp1");

        //fill_unfold_matrixs();
        //fill_unfold_fake_hists();


        /*
        double detector_dimass = (*p.lepton0 + *p.lepton1).M();
        vector<double> deltaR_bins;
        for (int i=0; i <= 500; i++){
            deltaR_bins.push_back(i*0.02);
        }
        // loop over photons
        // p.lepton0 p.lepton1
        double deltaR0 = gen_isr_l0.DeltaR(*p.lepton0) < gen_isr_l0.DeltaR(*p.lepton1)? gen_isr_l0.DeltaR(*p.lepton0) : gen_isr_l0.DeltaR(*p.lepton1);
        double deltaR1 = gen_isr_l1.DeltaR(*p.lepton0) < gen_isr_l1.DeltaR(*p.lepton1)? gen_isr_l1.DeltaR(*p.lepton0) : gen_isr_l1.DeltaR(*p.lepton1);
        FillHist(p.prefix+"deltaR_btw_detector_preFSR_lepton_detector_dimass", detector_dimass, deltaR0, p.weightmap[""],
                 bins.at("dimass").at("coarse").size()-1, bins.at("dimass").at("coarse").data(),
                 deltaR_bins.size()-1, deltaR_bins.data());
        FillHist(p.prefix+"deltaR_btw_detector_preFSR_lepton_detector_dimass", detector_dimass, deltaR1, p.weightmap[""],
                 bins.at("dimass").at("coarse").size()-1, bins.at("dimass").at("coarse").data(),
                 deltaR_bins.size()-1, deltaR_bins.data());

        deltaR0 = dressed_p1_l0.DeltaR(*p.lepton0) < dressed_p1_l0.DeltaR(*p.lepton1)? dressed_p1_l0.DeltaR(*p.lepton0) : dressed_p1_l0.DeltaR(*p.lepton1);
        deltaR1 = dressed_p1_l1.DeltaR(*p.lepton0) < dressed_p1_l1.DeltaR(*p.lepton1)? dressed_p1_l1.DeltaR(*p.lepton0) : dressed_p1_l1.DeltaR(*p.lepton1);
        FillHist(p.prefix+"deltaR_btw_detector_preFSR_p1_lepton_detector_dimass", detector_dimass, deltaR0, p.weightmap[""],
                 bins.at("dimass").at("coarse").size()-1, bins.at("dimass").at("coarse").data(),
                 deltaR_bins.size()-1, deltaR_bins.data());
        FillHist(p.prefix+"deltaR_btw_detector_preFSR_p1_lepton_detector_dimass", detector_dimass, deltaR1, p.weightmap[""],
                 bins.at("dimass").at("coarse").size()-1, bins.at("dimass").at("coarse").data(),
                 deltaR_bins.size()-1, deltaR_bins.data());
        */
    }
    fill_unfold_reco_hists();
    // lepton ID variable
    // lepton kinematics
    //FillLeptonHists(p);
}

void ISRAnalyzer::FillLeptonHists(Parameter& p){
    // array of dimass bins
    // array of ID variables
    const int n_mass_bin = 5;
    double dipt_cut = 100.0;
    double mass_edges[n_mass_bin + 1] = {55, 64, 81, 101, 200, 1000};
    double detector_dimass = (*p.lepton0 + *p.lepton1).M();
    double detector_dipt = (*p.lepton0 + *p.lepton1).Pt();

    if (detector_dipt < dipt_cut){
        for (int i = 0; i < n_mass_bin; i++){
            if (mass_edges[i] < detector_dimass && detector_dimass <= mass_edges[i + 1]) {
                // lepton kinematics
                TString hist_postfix = Form("_dipt_%.1f_dimass_%.1fto%.1f", dipt_cut, mass_edges[i], mass_edges[i+1]);
                FillHist(p.prefix + p.hprefix + "lep0_pt" + hist_postfix, (*p.lepton0).Pt(), p.weightmap[""], 500, 0., 500.);
                FillHist(p.prefix + p.hprefix + "lep1_pt" + hist_postfix, (*p.lepton1).Pt(), p.weightmap[""], 500, 0., 500.);
                FillHist(p.prefix + p.hprefix + "lep_pt" + hist_postfix, (*p.lepton0).Pt(), p.weightmap[""], 500, 0., 500.);
                FillHist(p.prefix + p.hprefix + "lep_pt" + hist_postfix, (*p.lepton1).Pt(), p.weightmap[""], 500, 0., 500.);

                FillHist(p.prefix + p.hprefix + "lep0_eta" + hist_postfix, (*p.lepton0).Eta(), p.weightmap[""], 120, -3.0, 3.0);
                FillHist(p.prefix + p.hprefix + "lep1_eta" + hist_postfix, (*p.lepton1).Eta(), p.weightmap[""], 120, -3.0, 3.0);
                FillHist(p.prefix + p.hprefix + "lep_eta" + hist_postfix, (*p.lepton0).Eta(), p.weightmap[""], 120, -3.0, 3.0);
                FillHist(p.prefix + p.hprefix + "lep_eta" + hist_postfix, (*p.lepton1).Eta(), p.weightmap[""], 120, -3.0, 3.0);

                FillHist(p.prefix + p.hprefix + "lep0_phi" + hist_postfix, (*p.lepton0).Phi(), p.weightmap[""], 160, -4.0, 4.0);
                FillHist(p.prefix + p.hprefix + "lep1_phi" + hist_postfix, (*p.lepton1).Phi(), p.weightmap[""], 160, -4.0, 4.0);
                FillHist(p.prefix + p.hprefix + "lep_phi" + hist_postfix, (*p.lepton0).Phi(), p.weightmap[""], 160, -4.0, 4.0);
                FillHist(p.prefix + p.hprefix + "lep_phi" + hist_postfix, (*p.lepton1).Phi(), p.weightmap[""], 160, -4.0, 4.0);

                FillHist(p.prefix + p.hprefix + "lep_phi" + add_postfix_BE_EE(p, 0, true) + hist_postfix, (*p.lepton0).Phi(), p.weightmap[""], 160, -4.0, 4.0);
                FillHist(p.prefix + p.hprefix + "lep_phi" + add_postfix_BE_EE(p, 1, true) + hist_postfix, (*p.lepton1).Phi(), p.weightmap[""], 160, -4.0, 4.0);

                if (p.channel == "ee") {
                    FillHist(p.prefix + p.hprefix + "lep_reliso" + add_postfix_BE_EE(p, 0) + hist_postfix, (*p.lepton0).RelIso(), p.weightmap[""], 200, 0.0, 0.2);
                    FillHist(p.prefix + p.hprefix + "lep_reliso" + add_postfix_BE_EE(p, 1) + hist_postfix, (*p.lepton1).RelIso(), p.weightmap[""], 200, 0.0, 0.2);


                    FillHist(p.prefix + p.hprefix + "lep_sigmaIetaIeta" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).Full5x5_sigmaIetaIeta(), p.weightmap[""], 100, 0., 0.05);
                    FillHist(p.prefix + p.hprefix + "lep_sigmaIetaIeta" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).Full5x5_sigmaIetaIeta(), p.weightmap[""], 100, 0., 0.05);

                    FillHist(p.prefix + p.hprefix + "lep_dEtaSeed" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).dEtaSeed(), p.weightmap[""], 200, -0.02, 0.02);
                    FillHist(p.prefix + p.hprefix + "lep_dEtaSeed" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).dEtaSeed(), p.weightmap[""], 200, -0.02, 0.02);

                    FillHist(p.prefix + p.hprefix + "lep_dPhiIn" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).dPhiIn(), p.weightmap[""], 200, -0.2, 0.2);
                    FillHist(p.prefix + p.hprefix + "lep_dPhiIn" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).dPhiIn(), p.weightmap[""], 200, -0.2, 0.2);

                    FillHist(p.prefix + p.hprefix + "lep_HoverE" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).HoverE(), p.weightmap[""], 200, 0., 0.2);
                    FillHist(p.prefix + p.hprefix + "lep_HoverE" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).HoverE(), p.weightmap[""], 200, 0., 0.2);

                    FillHist(p.prefix + p.hprefix + "lep_InvEminusInvP" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).InvEminusInvP(), p.weightmap[""], 200, 0., 0.2);
                    FillHist(p.prefix + p.hprefix + "lep_InvEminusInvP" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).InvEminusInvP(), p.weightmap[""], 200, 0., 0.2);

                    FillHist(p.prefix + p.hprefix + "lep_e2x5OverE5x5" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).e2x5OverE5x5(), p.weightmap[""], 200, 0., 2.0);
                    FillHist(p.prefix + p.hprefix + "lep_e2x5OverE5x5" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).e2x5OverE5x5(), p.weightmap[""], 200, 0., 2.0);

                    FillHist(p.prefix + p.hprefix + "lep_e1x5OverE5x5" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).e1x5OverE5x5(), p.weightmap[""], 200, 0., 2.0);
                    FillHist(p.prefix + p.hprefix + "lep_e1x5OverE5x5" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).e1x5OverE5x5(), p.weightmap[""], 200, 0., 2.0);

                    FillHist(p.prefix + p.hprefix + "lep_TrkIso" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).TrkIso(), p.weightmap[""], 100, 0., 100.0);
                    FillHist(p.prefix + p.hprefix + "lep_TrkIso" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).TrkIso(), p.weightmap[""], 100, 0., 100.0);

                    FillHist(p.prefix + p.hprefix + "lep_dr03EcalRecHitSumEt" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).dr03EcalRecHitSumEt(), p.weightmap[""], 100, 0., 100.0);
                    FillHist(p.prefix + p.hprefix + "lep_dr03EcalRecHitSumEt" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).dr03EcalRecHitSumEt(), p.weightmap[""], 100, 0., 100.0);

                    FillHist(p.prefix + p.hprefix + "lep_dr03HcalDepth1TowerSumEt" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).dr03HcalDepth1TowerSumEt(), p.weightmap[""], 100, 0., 100.0);
                    FillHist(p.prefix + p.hprefix + "lep_dr03HcalDepth1TowerSumEt" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).dr03HcalDepth1TowerSumEt(), p.weightmap[""], 100, 0., 100.0);

                    FillHist(p.prefix + p.hprefix + "lep_dr03HcalTowerSumEt" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).dr03HcalTowerSumEt(), p.weightmap[""], 100, 0., 100.0);
                    FillHist(p.prefix + p.hprefix + "lep_dr03HcalTowerSumEt" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).dr03HcalTowerSumEt(), p.weightmap[""], 100, 0., 100.0);

                    FillHist(p.prefix + p.hprefix + "lep_dr03TkSumPt" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).dr03TkSumPt(), p.weightmap[""], 100, 0., 100.0);
                    FillHist(p.prefix + p.hprefix + "lep_dr03TkSumPt" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).dr03TkSumPt(), p.weightmap[""], 100, 0., 100.0);

                    FillHist(p.prefix + p.hprefix + "lep_ecalPFClusterIso" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).ecalPFClusterIso(), p.weightmap[""], 100, 0., 100.0);
                    FillHist(p.prefix + p.hprefix + "lep_ecalPFClusterIso" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).ecalPFClusterIso(), p.weightmap[""], 100, 0., 100.0);

                    FillHist(p.prefix + p.hprefix + "lep_hcalPFClusterIso" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).hcalPFClusterIso(), p.weightmap[""], 100, 0., 100.0);
                    FillHist(p.prefix + p.hprefix + "lep_hcalPFClusterIso" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).hcalPFClusterIso(), p.weightmap[""], 100, 0., 100.0);

                    FillHist(p.prefix + p.hprefix + "lep_isEcalDriven" + add_postfix_BE_EE(p, 0) + hist_postfix, (p.electrons.at(0)).isEcalDriven(), p.weightmap[""], 2, 0., 2.0);
                    FillHist(p.prefix + p.hprefix + "lep_isEcalDriven" + add_postfix_BE_EE(p, 1) + hist_postfix, (p.electrons.at(1)).isEcalDriven(), p.weightmap[""], 2, 0., 2.0);
                }
            } // mass edges
        }
    } // dipt cut
}

TString ISRAnalyzer::add_postfix_BE_EE(Parameter& p, int lepton_index, bool add_sign){
    double eta;
    double eta_split;
    TString lepton_name;
    if (p.channel == "ee") {
        lepton_name = "E";
        eta_split = 1.479;
        eta = p.electrons.at(lepton_index).scEta();
    } else {
        lepton_name = "M";
        eta_split = 0.9;
        eta = p.muons.at(lepton_index).Eta();
    }

    TString postfix;
    if (fabs(eta) <= eta_split) {
        postfix = "_B" + lepton_name;
    }
    else {
        postfix = "_E" + lepton_name;
    }

    if (add_sign) {
        if (eta >= 0)
            postfix += "+";
        else
            postfix += "-";
    }
    return postfix;
}

void ISRAnalyzer::get_DY_gen_particles(){
    if(!IsDYSample){
        cout <<"[ISRAnalyzer::get_DY_gen_particles] this is for DY event"<<endl;
        exit(EXIT_FAILURE);
    }

    int ngen=gens.size();
    // set parton0, parton1
    // set leptons, photons
    // collect status 1 leptons and photons
    for(int i = 0; i < ngen; i++){
        if(!gens.at(i).isPrompt()) continue;

        int pid = gens.at(i).PID();
        int mother_index  = gens.at(i).MotherIndex();
        if(mother_index == -1) continue;
        //int mother_pid = gens.at(mother_index).PID();
        /*
        if(gens.at(i).Status() == 1 || gens.at(i).isPromptFinalState() == true)
        {
            std::cout << i << "\033[1;31m pid: " << gens.at(i).PID()
                           << " mother index: " << gens.at(i).MotherIndex()
                           << " status: " << gens.at(i).Status()
                           << " isPrompt: " << gens.at(i).isPrompt()
                           << " isHardProcess: " << gens.at(i).isHardProcess()
                           << " isPromptFinalState: " << gens.at(i).isPromptFinalState()
                           << " Pt: " << gens.at(i).Pt()
                           << "\033[0m" << std::endl;
        }
        else
        {
            std::cout << i << " pid: " << gens.at(i).PID() << " mother index: " << gens.at(i).MotherIndex()
                       << " status: " << gens.at(i).Status() << " isPrompt: " << gens.at(i).isPrompt()
                       << " isHardProcess: " << gens.at(i).isHardProcess()
                       << " Pt: " << gens.at(i).Pt()
                       << std::endl;
        }
        */
        if(gens.at(i).isHardProcess()){
            if(abs(pid) < TOP || pid == PHOTON){
                if(gen_p0_isr.IsEmpty()){
                    gen_p0_isr = gens[i];
                }
                else{
                    if(gen_p1_isr.IsEmpty())
                        gen_p1_isr = gens[i];
                }
            }
        }// hard process

        if(gens.at(i).Status() == 1){
            //if(mother_pid != PROTON){ // stable lepton from proton, gamma gamma?
            switch(abs(pid)) {
                case ELECTRON :
                case MUON :
                    gen_leptons.push_back(&gens[i]);
                    break;
                case PHOTON :
                    gen_photons.push_back(&gens[i]);
                    break;
                default :
                    // do nothing
                    break;
            }
            //}
        }// status 1
    }// gen loop
    get_DY_bare_lepton_pair();
}

//
void ISRAnalyzer::get_DY_bare_lepton_pair(bool verbose){


    int nlepton=gen_leptons.size();

    for(int i = 0; i < nlepton; i++){

        vector<int> lepton0_history;
        save_gen_history(gens, *gen_leptons[i], lepton0_history);

        for(int j = i + 1; j < nlepton; j++){

            if(!(gen_leptons[i]->PID()+gen_leptons[j]->PID() == 0)) continue;

            vector<int> lepton1_history;
            save_gen_history(gens, *gen_leptons[j], lepton1_history);

            vector<int> history_intersection;
            std::set_intersection(lepton0_history.begin(), lepton0_history.end(),
                                  lepton1_history.begin(), lepton1_history.end(),
                                  back_inserter(history_intersection));

            if(history_intersection.size() == 0) continue;
            std::sort(history_intersection.begin(), history_intersection.end(), greater <>());

            // check if lepton vertex included ex) e^- to e^-e^-e^+
            bool letpon_vertex_included = false;
            for(auto vertex : history_intersection){
                int pid_vertex = gens.at(vertex).PID();
                if(abs(pid_vertex) == ELECTRON || abs(pid_vertex) == MUON){
                    letpon_vertex_included = true;
                    break;
                }
            }
            if(letpon_vertex_included) continue;
            //DY_index = history_intersection.at(0);
            if(verbose){
                cout << "i,j " << i << "," << j << endl;
                cout << "(lepton0 + lepton1).M() : " << (gen_l0_bare + gen_l1_bare).M() << endl;
                cout << "(*gen_leptons[i] + *gen_leptons[j]).M() : " << (*gen_leptons[i] + *gen_leptons[j]).M() << endl;
                cout << gen_leptons[i]->Pt() << " " << gen_leptons[j]->Pt() << endl;
            }

            // set DY lepton pair
            if((*gen_leptons[i] + *gen_leptons[j]).M() > (gen_l0_bare + gen_l1_bare).M()){ // note that there could be negative invariant mass event
                DY_index = history_intersection.at(0); // set DY index when gen_leptons selected
                if(verbose)
                    cout << "set lepton!" << endl;

                if(gen_leptons[i]->Pt()>gen_leptons[j]->Pt()){
                    gen_l0_bare = *gen_leptons[i];
                    gen_l1_bare = *gen_leptons[j];
                    gen_l0_bare.SetIndexPIDStatus(gen_leptons[i]->Index(), gen_leptons[i]->PID(), gen_leptons[i]->Status());
                    gen_l1_bare.SetIndexPIDStatus(gen_leptons[j]->Index(), gen_leptons[j]->PID(), gen_leptons[j]->Status());
                }
                else{
                    gen_l0_bare = *gen_leptons[j];
                    gen_l1_bare = *gen_leptons[i];
                    gen_l0_bare.SetIndexPIDStatus(gen_leptons[j]->Index(), gen_leptons[j]->PID(), gen_leptons[j]->Status());
                    gen_l1_bare.SetIndexPIDStatus(gen_leptons[i]->Index(), gen_leptons[i]->PID(), gen_leptons[i]->Status());
                }
            }
        }// leptons(j) loop
    }// leptons(i) loop
}

void ISRAnalyzer::get_DY_dressed_lepton_pair(Gen& lepton0, Gen& lepton1, const DressedMode mode, const double dR){

    if(DY_index == -1) // DY leptons not selected
        return;

    vector<int> DY_history;
    save_gen_history(gens, gen_l0_bare, DY_history, DY_index);
    save_gen_history(gens, gen_l1_bare, DY_history, DY_index);
    DY_history.push_back(DY_index);

    lepton0 = gen_l0_bare;
    lepton1 = gen_l1_bare;

    // status 1 photons
    for(const auto photon : gen_photons){
        if(mode == DressedMode::MotherMatch || mode == DressedMode::MotherDRMatch){
            vector<int> history_photon = TrackGenSelfHistory(*photon, gens);
            auto it = find(DY_history.begin(), DY_history.end(), history_photon.at(1));
            // photon's mother not exist in DY history
            if(it == DY_history.end()){
                continue;
            }
        }
        if(mode == DressedMode::DRMatch || mode == DressedMode::MotherDRMatch){
            if(gen_l0_bare.DeltaR(*photon) < gen_l1_bare.DeltaR(*photon)){
                if(gen_l0_bare.DeltaR(*photon) > dR) continue;
            }
            else{
                if(gen_l1_bare.DeltaR(*photon) > dR) continue;
            }
        }
        //added_photons.push_back(photon);

        // add gamma to the cloest lepton
        if(gen_l0_bare.DeltaR(*photon) < gen_l1_bare.DeltaR(*photon)){
            lepton0 +=  *photon;
        }
        else{
            lepton1 +=  *photon;
        }
    }// loops photon
}

// save index starting from the current index of particle to the final index of mother particle
void ISRAnalyzer::save_gen_history(const vector<Gen>& gens, const Gen& lepton, vector<int>& index_vector, const int index_limit){

    //int pid = lepton.PID();
    int index = lepton.Index();
    int mother_index  = lepton.MotherIndex();
    index_vector.push_back(index);

    while(mother_index > index_limit){ //
        index_vector.push_back(mother_index);
        index = mother_index;
        mother_index = gens.at(index).MotherIndex();
    }
    sort(index_vector.begin(), index_vector.end());
}

void ISRAnalyzer::print_gen_particles(const vector<Gen>& gens){

    for(unsigned int i = 0; i < gens.size(); i++){
        if(gens.at(i).Status() == 1 || gens.at(i).isPromptFinalState() == true)
        {
            std::cout << i << "\033[1;31m pid: " << gens.at(i).PID()
            << " mother index: " << gens.at(i).MotherIndex()
            << " status: " << gens.at(i).Status()
            << " isPrompt: " << gens.at(i).isPrompt()
            << " isHardProcess: " << gens.at(i).isHardProcess()
            << " isPromptFinalState: " << gens.at(i).isPromptFinalState()
            << " Pt: " << gens.at(i).Pt()
            << "\033[0m" << std::endl;
        }
        else
        {
            std::cout << i << " pid: " << gens.at(i).PID() << " mother index: " << gens.at(i).MotherIndex()
            << " status: " << gens.at(i).Status() << " isPrompt: " << gens.at(i).isPrompt()
            << " isHardProcess: " << gens.at(i).isHardProcess()
            << " Pt: " << gens.at(i).Pt()
            << std::endl;
        }
    }
}

ISRAnalyzer::ISRAnalyzer(){
    job_number=-1;
    // 2D
    ISRUnfoldBin* dipt_2d_fine = create_2d_unfold_bin("dipt", "fine", false, true, "dimass", "window_v1", true, true);
    ISRUnfoldBin* dipt_2d_coarse = create_2d_unfold_bin("dipt", "coarse", false, true, "dimass", "window_v1", true, true);
    ISRUnfoldBin* dimass_2d_fine = create_2d_unfold_bin("dimass", "fine", true, true, "dipt", "window_v1", false, true);
    ISRUnfoldBin* dimass_2d_coarse = create_2d_unfold_bin("dimass", "coarse", true, true, "dipt", "window_v1", false, true);

    //create_2d_unfold_set(dipt_2d_fine, dipt_2d_fine, 0, 1500, 50, 1500);
    create_2d_unfold_set(dipt_2d_fine, dipt_2d_coarse, 0, 1500, 50, 1500);
    //create_2d_unfold_set(dimass_2d_fine, dimass_2d_fine, 0, 1500, 50, 1500);
    //create_2d_unfold_set(dimass_2d_fine, dimass_2d_coarse, 0, 1500, 50, 1500, true, false);

    // 1D
    vector<double> mass_window = {55, 1500};  // DO NOT ALLOW UNDERFLOW event for 1D mass unfolding
    ISRUnfoldBin* dimass_1d_fine = create_1d_unfold_bin("dimass", "fine", false, true);  // for 1D, underflow is the 0th bin!!
    ISRUnfoldBin* dimass_1d_coarse = create_1d_unfold_bin("dimass", "coarse", false, true);
    create_1d_unfold_set(dimass_1d_fine, dimass_1d_coarse, 0, 100.0, mass_window);
    create_1d_unfold_set(dimass_1d_coarse, dimass_1d_coarse, 0, 100.0, mass_window, false, true);

    double dipt_cut = 150;
    ISRUnfoldBin* dipt_1d_fine = create_1d_unfold_bin("dipt", "fine", false, true);
    ISRUnfoldBin* dipt_1d_coarse = create_1d_unfold_bin("dipt", "coarse", false, true);

    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 55, 64);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 64, 81);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 81, 101);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 101, 200);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 200, 1000);

/*
    ISRUnfoldBin* dipt_1d_v3 = create_1d_unfold_bin("dipt", "coarse_v3");

    ISRUnfoldBin* dipt_1d_ext_v1 = create_1d_unfold_bin("dipt", "fine_extended");
    ISRUnfoldBin* dipt_1d_ext_v2 = create_1d_unfold_bin("dipt", "coarse_extended");
    ISRUnfoldBin* dipt_1d_ext_v3 = create_1d_unfold_bin("dipt", "coarse_v3_extended");

    vector<double> mass_windows = {55, 64, 81, 101, 200, 1000};

    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 55, 68);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 55, 81);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 64, 81);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 72, 91);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 81, 101);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 91, 110);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 96, 106);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 96, 126);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 101, 150);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 101, 200);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 106, 220);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 110, 243);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 115, 273);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 120, 320);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 126, 380);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 133, 440);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 141, 510);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 150, 600);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 160, 700);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 171, 830);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 185, 1000);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_coarse, 0, dipt_cut, 200, 1000);

    create_1d_unfold_set(dipt_1d_fine, dipt_1d_v3, 0, dipt_cut, 126, 380, true, false);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_v3, 0, dipt_cut, 133, 440, true, false);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_v3, 0, dipt_cut, 141, 510, true, false);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_v3, 0, dipt_cut, 150, 600, true, false);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_v3, 0, dipt_cut, 160, 700, true, false);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_v3, 0, dipt_cut, 171, 830, true, false);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_v3, 0, dipt_cut, 185, 1000, true, false);
    create_1d_unfold_set(dipt_1d_fine, dipt_1d_v3, 0, dipt_cut, 200, 1000, true, false);

    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 55, 68, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 55, 81, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 72, 91, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 91, 110, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 96, 106, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 96, 126, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 101, 150, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 106, 220, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 110, 243, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 115, 273, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 120, 320, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 126, 380, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 133, 440, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 141, 510, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 150, 600, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 160, 700, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 171, 830, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_coarse, 0, dipt_cut, 185, 1000, false, true);

    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_v3, 0, dipt_cut, 126, 380, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_v3, 0, dipt_cut, 133, 440, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_v3, 0, dipt_cut, 141, 510, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_v3, 0, dipt_cut, 150, 600, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_v3, 0, dipt_cut, 160, 700, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_v3, 0, dipt_cut, 171, 830, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_v3, 0, dipt_cut, 185, 1000, false, true);
    create_1d_unfold_set(dipt_1d_coarse, dipt_1d_v3, 0, dipt_cut, 200, 1000, false, true);

    // extended dipt cut
    dipt_cut = 500;
    // fine coarse
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 55, 64);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 55, 68);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 55, 81);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 64, 81);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 72, 91);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 81, 101);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 91, 110);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 96, 106);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 96, 126);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 101, 150);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 101, 200);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 106, 220);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 110, 243);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 115, 273);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 120, 320);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 126, 380);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 133, 440);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 141, 510);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 150, 600);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 160, 700);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 171, 830);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 185, 1000);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v2, 0, dipt_cut, 200, 1000);

    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 55, 64, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 55, 68, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 55, 81, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 64, 81, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 72, 91, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 126, 380, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 133, 440, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 141, 510, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 150, 600, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 160, 700, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 171, 830, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 185, 1000, true, false);
    create_1d_unfold_set(dipt_1d_ext_v1, dipt_1d_ext_v3, 0, dipt_cut, 200, 1000, true, false);

    create_1d_unfold_set(dimass_1d_fine, dimass_1d_coarse, 0, 150.0, mass_window);
    create_1d_unfold_set(dimass_1d_fine, dimass_1d_coarse, 0, 300.0, mass_window);
    create_1d_unfold_set(dimass_1d_fine, dimass_1d_coarse, 0, 1000.0, mass_window);

    create_1d_unfold_set(dimass_1d_coarse, dimass_1d_coarse, 0, 100.0, mass_window, false, true);
    create_1d_unfold_set(dimass_1d_coarse, dimass_1d_coarse, 0, 150.0, mass_window, false, true);
    create_1d_unfold_set(dimass_1d_coarse, dimass_1d_coarse, 0, 300.0, mass_window, false, true);
    create_1d_unfold_set(dimass_1d_coarse, dimass_1d_coarse, 0, 1000.0, mass_window, false, true);
    */
}
ISRAnalyzer::~ISRAnalyzer(){
    //DeleteCosThetaWeight();
}
