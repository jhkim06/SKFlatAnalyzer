#include "ISRAnalyzer.h"

void ISRAnalyzer::initializeAnalyzer(){
    
    SMPAnalyzerCore::initializeAnalyzer(); //setup zpt roc z0
    ISRUnfold::initializeISRUnfold(job_number);
    
    vector<JetTagging::Parameters> jtps={JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb)};
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
}

void ISRAnalyzer::executeEvent(){
    // FIXME some events of DYJets has nan PDF weights. I don't know why...
    if(MCSample=="DYJets"&&!isnormal(weight_Scale->at(0))) return;
    // GEN level
    executeEventGen();
    
    // RECO level
    if(!IsDATA||DataStream.Contains("DoubleMuon")){
        executeEventWithParameter(MakeParameter("mm"));
    }
    if(!IsDATA||DataStream.Contains("DoubleEG")||DataStream.Contains("EGamma")){
        executeEventWithParameter(MakeParameter("ee"));
    }
}

SMPAnalyzerCore::Parameter ISRAnalyzer::MakeParameter(TString key){
    
    Parameter p = SMPAnalyzerCore::MakeParameter(key);
    if (key == "ee"){
        p.SetElectrons(SMPGetElectrons("passMediumID",0.0,2.5));
    }
    if (key == "mm"){
        p.SetMuons(MuonMomentumCorrection(SMPGetMuons("POGTightWithTightIso",0.0,2.4),0,0));  
    }
    
    p.weightbit = 0;
    
    if(IsNominalRun) p.weightbit|=NominalWeight;
    if(HasFlag("SYS")&&!IsDATA&&(p.channel=="ee"||p.channel=="mm")) p.weightbit|=SystematicWeight|EfficiencyWeight;
    if(HasFlag("PDFSYS")&&!IsDATA&&(p.channel=="ee"||p.channel=="mm")) p.weightbit|=PDFWeight;
    
    if(HasFlag("nbjet")) p.prefix+="nbjet/";
    else if(HasFlag("0bjet")) p.prefix+="0bjet/";
    if(HasFlag("highmet")) p.prefix+="highmet/";
    
    return p;
}

// called in SMPAnalyzerCore::executeEventWithParameter()
bool ISRAnalyzer::PassSelection(Parameter& p){
    
    int n_bjet=0;
    std::vector<Jet> jets=GetJets("tightLepVeto",30,2.4);
    std::sort(jets.begin(),jets.end(),PtComparing);
    JetTagging::Parameters jtp = JetTagging::Parameters(JetTagging::DeepJet,JetTagging::Tight,JetTagging::incl,JetTagging::comb);
    for(const auto& jet:jets)
        if(jet.GetTaggerResult(jtp.j_Tagger) > mcCorr->GetJetTaggingCutValue(jtp.j_Tagger, jtp.j_WP))
            n_bjet++;
    
    p.intmap["nbjet"]=n_bjet;
    p.doublemap["btagSF"]=mcCorr->GetBTaggingReweight_1a(jets,jtp);
    p.doublemap["btagSF_hup"]=mcCorr->GetBTaggingReweight_1a(jets,jtp,"SystUpHTag");
    p.doublemap["btagSF_hdown"]=mcCorr->GetBTaggingReweight_1a(jets,jtp,"SystDownHTag");
    p.doublemap["btagSF_lup"]=mcCorr->GetBTaggingReweight_1a(jets,jtp,"SystUpLTag");
    p.doublemap["btagSF_ldown"]=mcCorr->GetBTaggingReweight_1a(jets,jtp,"SystDownLTag");
    
    if(n_bjet > 0) return false; // reguire bjet veto
    if(p.prefix.Contains("nbjet")&&!n_bjet) return false;
    if(p.prefix.Contains("0bjet")&&n_bjet) return false;

    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCut",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight);
    if(IsNominalRun) FillCutflow(p.prefix+p.hprefix+"cutflow","BJetCutSF",p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.doublemap["btagSF"]);
    
    if(!SMPAnalyzerCore::PassSelection(p)) return false;
    return true;
    
}

void ISRAnalyzer::executeEventGen(){
    
    if(IsDYSample){
        if(abs(lhe_l0.ID())!=15&&abs(lhe_l1.ID())!=15){
            Parameter p;
            if( (abs(lhe_l0.ID())==11&&abs(lhe_l1.ID())==11) || (!lhes.size()&&abs(gen_l0.PID())==11&&abs(gen_l1.PID())==11) ){
                p=MakeParameter("ee");
            } else if( (abs(lhe_l0.ID())==13&&abs(lhe_l1.ID())==13) || (!lhes.size()&&abs(gen_l0.PID())==13&&abs(gen_l1.PID())==13) ){
                p=MakeParameter("mm");
            } else {
                exit(EXIT_FAILURE);
            }
            
            const vector<Gen> gens=GetGens();
            Gen gen_isr_parton0, gen_isr_parton1;
            Gen gen_isr_l0, gen_isr_l1; // dressed gen particles
            vector<const Gen*> added_photons;
            
            int DY_index = get_DY_gen_particles(gens, gen_isr_parton0, gen_isr_parton1, gen_isr_l0, gen_isr_l1, PreFSR, added_photons);
            //SMPAnalyzerCore::GetAFBGenParticles(gens, gen_isr_parton0, gen_isr_parton1, gen_isr_l0, gen_isr_l1, 3);
            TLorentzVector dilepton = gen_isr_l0 + gen_isr_l1;
            double dimass=dilepton.M();
            double dirap=dilepton.Rapidity();
            double dipt=dilepton.Pt();
            // apply only dipt and dimass cut
            
            map<TString,double> map_weight;
            map_weight[""]=p.w.lumiweight; // TODO check which gen lepton used to get weight
            
            if (dipt < 1500 && dimass > 53 && dimass < 1500) {
                fill_unfold_hists(p, (Particle*)&gen_isr_l0, (Particle*)&gen_isr_l1,
                                  map_weight, TUnfoldBin::unfolded_bin, "gen_acceptance");
            }
        }
    }
}

void ISRAnalyzer::executeEventWithParameter(Parameter& p){
    // draw cuflow according to the lepton selection defined in p
    SMPAnalyzerCore::executeEventWithParameter(p);
}

void ISRAnalyzer::EvalWeights(Parameter& p){
    
    if(p.weightbit&NominalWeight){
        // make function to get weight
        p.weightmap[""]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
    }
    if(p.weightbit&SystematicWeight){
        if(!IsDATA){
            p.weightmap["_zptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.zptweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
            p.weightmap["_noPUweight"]=p.w.lumiweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]; //need for AN
            p.weightmap["_PUweight_up"]=p.w.lumiweight*p.w.PUweight_up*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
            p.weightmap["_PUweight_down"]=p.w.lumiweight*p.w.PUweight_down*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
            
            p.weightmap["_noprefireweight"]=p.w.lumiweight*p.w.PUweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
            p.weightmap["_prefireweight_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_up*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
            p.weightmap["_prefireweight_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight_down*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
            
            //p.weightmap["_nozptweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
            p.weightmap["_noz0weight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
            p.weightmap["_noweakweight"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
            
            p.weightmap["_nobtagSF"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF;
            p.weightmap["_btagSF_hup"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF_hup"];
            p.weightmap["_btagSF_hdown"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF_hdown"];
            p.weightmap["_btagSF_lup"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF_lup"];
            p.weightmap["_btagSF_ldown"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF_ldown"];
            
            for(int j=0,nj=fEff->nreplica;j<nj;j++){
                double electronRECOSF=p.w.electronRECOSF_sys.size() ? p.w.electronRECOSF_sys[0][j] : 1.;
                double electronIDSF=p.w.electronIDSF_sys.size() ? p.w.electronIDSF_sys[0][j] : 1.;
                double muonIDSF=p.w.muonIDSF_sys.size() ? p.w.muonIDSF_sys[0][j] : 1.;
                double triggerSF=p.w.triggerSF_sys.size() ? p.w.triggerSF_sys[0][j] : 1.;
                p.weightmap[Form("_efficiencySF_stat%d",j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*electronRECOSF*electronIDSF*muonIDSF*p.w.muonISOSF*triggerSF*p.doublemap["btagSF"];
            }
            
            for(int i=1,ni=p.w.electronRECOSF_sys.size();i<ni;i++){
                for(int j=0,nj=p.w.electronRECOSF_sys[i].size();j<nj;j++){
                    p.weightmap[Form("_electronRECOSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF_sys[i][j]*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
                }
            }
            
            for(int i=1,ni=p.w.electronIDSF_sys.size();i<ni;i++){
                for(int j=0,nj=p.w.electronIDSF_sys[i].size();j<nj;j++){
                    p.weightmap[Form("_electronIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF_sys[i][j]*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
                }
            }
            
            for(int i=1,ni=p.w.muonIDSF_sys.size();i<ni;i++){
                for(int j=0,nj=p.w.muonIDSF_sys[i].size();j<nj;j++){
                    p.weightmap[Form("_muonIDSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF_sys[i][j]*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"];
                }
            }
            
            for(int i=1,ni=p.w.triggerSF_sys.size();i<ni;i++){
                for(int j=0,nj=p.w.triggerSF_sys[i].size();j<nj;j++){
                    p.weightmap[Form("_triggerSF_s%d_m%d",i,j)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF_sys[i][j]*p.doublemap["btagSF"];
                }
            }
        }
    }
    if(p.weightbit&PDFWeight){
        for(unsigned int i=0;i<weight_Scale->size();i++){
            p.weightmap[Form("_scalevariation%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_Scale->at(i);
        }
        for(unsigned int i=0;i<weight_PDF->size();i++){
            p.weightmap[Form("_pdf%d",i)]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_PDF->at(i);
        }
        if(weight_AlphaS->size()==2){
            p.weightmap["_alphaS_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_AlphaS->at(0);
            p.weightmap["_alphaS_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_AlphaS->at(1);
        }
        
        if(MCSample.Contains("MiNNLO")){
            p.weightmap["_sthw2_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_sthw2->at(0);
            p.weightmap["_sthw2_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_sthw2->at(2);
            p.weightmap["_largeptscales"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_largeptscales->at(0);
            p.weightmap["_q0_up"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_q0->at(0);
            p.weightmap["_q0_down"]=p.w.lumiweight*p.w.PUweight*p.w.prefireweight*p.w.z0weight*p.w.weakweight*p.w.electronRECOSF*p.w.electronIDSF*p.w.muonIDSF*p.w.muonISOSF*p.w.triggerSF*p.doublemap["btagSF"]*weight_q0->at(2);
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
    
    TLorentzVector dilepton=*p.lepton0+*p.lepton1;
    double dimass=dilepton.M();
    double dipt=dilepton.Pt();
    
    if (dipt < 1500 && dimass > 53 && dimass < 1500) {
        if (IsDYSample && p.hprefix!="tau_"){
            const vector<Gen> gens=GetGens();
            Gen gen_isr_parton0, gen_isr_parton1;
            Gen gen_isr_l0, gen_isr_l1;  // dressed gen particles
            Gen gen_bare_l0, gen_bare_l1;  // post FSR gen particles
            vector<const Gen*> added_photons;
            
            int DY_index = get_DY_gen_particles(gens, gen_isr_parton0, gen_isr_parton1, gen_isr_l0, gen_isr_l1, PreFSR, added_photons);
            DY_index = get_DY_gen_particles(gens, gen_isr_parton0, gen_isr_parton1, gen_bare_l0, gen_bare_l1);
            Parameter pgen=p;
            ResetRecoWeights(pgen);
            EvalWeights(pgen);

            // require the same kinematic cuts for the unfolded level as for the reco level
            if (pass_lepton_kinematic_selections(p, &gen_isr_l0, &gen_isr_l1)) {
                // unfolded level
                fill_unfold_hists(p, (Particle*)&gen_isr_l0, (Particle*)&gen_isr_l1,
                                  pgen.weightmap, TUnfoldBin::unfolded_bin, "gen_dRp1");
                // fill response matrix
                fill_unfold_response_matrixs(p, (Particle*)p.lepton0, (Particle*)p.lepton1, (Particle*)&gen_isr_l0, (Particle*)&gen_isr_l1,
                        p.weightmap, pgen.weightmap, "reco", "gen_dRp1");
            } else {
                // fake DY events
                fill_unfold_hists(p, (Particle*)p.lepton0, (Particle*)p.lepton1,
                                  p.weightmap, TUnfoldBin::folded_bin, "reco_fake");
            }
        }
        // reco level
        fill_unfold_hists(p, (Particle*)p.lepton0, (Particle*)p.lepton1,
                          p.weightmap, TUnfoldBin::folded_bin, "reco");
    }
}

bool ISRAnalyzer::pass_lepton_kinematic_selections(const Parameter& p, Particle* l0, Particle* l1){
    
    TLorentzVector dilepton=*l0+*l1;
    double dimass=dilepton.M();
    double dipt=dilepton.Pt();
    
    bool passed = false;
    if (( (*l0).Pt() > p.c.lepton0pt && (*l1).Pt( ) > p.c.lepton1pt)
        || ( (*l0).Pt() > p.c.lepton1pt && (*l1).Pt() > p.c.lepton0pt)) {
        
        double eta_cut = 2.5;
        if (p.channel.Contains("mm")) eta_cut = 2.4;
        
        if (fabs((*l0).Eta()) < eta_cut && fabs((*l1).Eta()) < eta_cut) {
            
            if (dipt < 1500 && dimass > 53 && dimass < 1500) {
                passed = true;
            }
        }
    }
    return passed;
}

int ISRAnalyzer::get_DY_gen_particles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& lepton0, Gen& lepton1, int mode, vector<const Gen*>& added_photons){
    
    int DY_index = -1;
    
    if(!IsDYSample){
        cout <<"[ISRAnalyzer::get_DY_gen_particles] this is for DY event"<<endl;
        exit(EXIT_FAILURE);
    }
    
    int ngen=gens.size();
    vector<const Gen*> leptons;
    vector<const Gen*> photons;
    
    // set parton0, parton1
    // set leptons, photons
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
                if(parton0.IsEmpty()){
                    parton0 = gens[i];
                }
                else{
                    if(parton1.IsEmpty())
                        parton1 = gens[i];
                }
            }
        }// hard process
        
        if(gens.at(i).Status() == 1){
            //if(mother_pid != PROTON){ // stable lepton from proton, gamma gamma?
            switch(abs(pid)) {
                case ELECTRON :
                case MUON :
                    leptons.push_back(&gens[i]);
                    break;
                case PHOTON :
                    photons.push_back(&gens[i]);
                    break;
                default :
                    // do nothing
                    break;
            }
            //}
        }// status 1
    }// gen loop
    
    // DY pair
    switch(mode) {
        case DressedDRp1 :
            DY_index = get_DY_dressed_lepton_pair(gens, leptons, photons, lepton0, lepton1, DressedMode::DRMatch, added_photons, 0.1);
            break;
        case DressedDRp4 :
            DY_index = get_DY_dressed_lepton_pair(gens, leptons, photons, lepton0, lepton1, DressedMode::DRMatch, added_photons, 0.4);
            break;
        case PreFSR :
            DY_index = get_DY_dressed_lepton_pair(gens, leptons, photons, lepton0, lepton1, DressedMode::MotherMatch, added_photons);
            break;
        default :
            break;
    }
    
    return DY_index;
}

int ISRAnalyzer::get_DY_gen_particles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& lepton0, Gen& lepton1){
    
    
    if(!IsDYSample){
        cout <<"[ISRAnalyzer::get_DY_gen_particles] this is for DY event"<<endl;
        exit(EXIT_FAILURE);
    }
    
    int ngen=gens.size();
    vector<const Gen*> leptons;
    vector<const Gen*> photons;
    
    for(int i = 0; i < ngen; i++){
        if(!gens.at(i).isPrompt()) continue;
        
        int pid = gens.at(i).PID();
        int mother_index  = gens.at(i).MotherIndex();
        if(mother_index == -1) continue;
        //int mother_pid = gens.at(mother_index).PID();
        
        if(gens.at(i).isHardProcess()){
            if(abs(pid) < TOP || pid == PHOTON){
                if(parton0.IsEmpty()){
                    parton0 = gens[i];
                }
                else{
                    if(parton1.IsEmpty())
                        parton1 = gens[i];
                }
            }
        }// hard process
        
        if(gens.at(i).Status() == 1){
            //if(mother_pid != PROTON){ // stable lepton from proton, gamma gamma?
            switch(abs(pid)) {
                case ELECTRON :
                case MUON :
                    leptons.push_back(&gens[i]);
                    break;
                case PHOTON :
                    photons.push_back(&gens[i]);
                    break;
                default :
                    // do nothing
                    break;
            }
            //}
        }// status 1
    }// gen loop
    
    int DY_index = get_DY_bare_lepton_pair(gens, leptons, lepton0, lepton1);
    return DY_index;
}

//
int ISRAnalyzer::get_DY_dressed_lepton_pair(const vector<Gen>& gens, const vector<const Gen*>& leptons, vector<const Gen*>& photons, Gen& lepton0, Gen& lepton1,
                                            const DressedMode mode, vector<const Gen*>& added_photons, const double dR){
    
    int DY_index = get_DY_bare_lepton_pair(gens, leptons, lepton0, lepton1);
    if(DY_index == -1) // DY leptons not selected
        return DY_index;
    
    vector<int> DY_history;
    save_gen_history(gens, lepton0, DY_history, DY_index);
    save_gen_history(gens, lepton1, DY_history, DY_index);
    DY_history.push_back(DY_index);
    
    Gen lepton0_temp = lepton0;
    Gen lepton1_temp = lepton1;
    
    // status 1 photons
    for(const auto photon : photons){
        if(mode == DressedMode::MotherMatch || mode == DressedMode::MotherDRMatch){
            vector<int> history_photon = TrackGenSelfHistory(*photon, gens);
            auto it = find(DY_history.begin(), DY_history.end(), history_photon.at(1));
            // photon's mother not exist in DY history
            if(it == DY_history.end()){
                continue;
            }
        }
        if(mode == DressedMode::DRMatch || mode == DressedMode::MotherDRMatch){
            if(lepton0_temp.DeltaR(*photon) < lepton1_temp.DeltaR(*photon)){
                if(lepton0_temp.DeltaR(*photon) > dR) continue;
            }
            else{
                if(lepton1_temp.DeltaR(*photon) > dR) continue;
            }
        }
        added_photons.push_back(photon);
        
        // add gamma to the cloest lepton
        if(lepton0_temp.DeltaR(*photon) < lepton1_temp.DeltaR(*photon)){
            lepton0 +=  *photon;
        }
        else{
            lepton1 +=  *photon;
        }
    }// loops photon
    return DY_index;
}


int ISRAnalyzer::get_DY_bare_lepton_pair(const vector<Gen>& gens, const vector<const Gen*>& leptons, Gen& lepton0, Gen& lepton1, bool verbose){
    
    int nlepton=leptons.size();
    int DY_index = -1;
    
    for(int i = 0; i < nlepton; i++){
        
        vector<int> lepton0_history;
        save_gen_history(gens, *leptons[i], lepton0_history);
        
        for(int j = i + 1; j < nlepton; j++){
            
            if(!(leptons[i]->PID()+leptons[j]->PID() == 0)) continue;
            
            vector<int> lepton1_history;
            save_gen_history(gens, *leptons[j], lepton1_history);
            
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
                cout << "(lepton0 + lepton1).M() : " << (lepton0 + lepton1).M() << endl;
                cout << "(*leptons[i] + *leptons[j]).M() : " << (*leptons[i] + *leptons[j]).M() << endl;
                cout << leptons[i]->Pt() << " " << leptons[j]->Pt() << endl;
            }
            
            // set DY lepton pair
            if((*leptons[i] + *leptons[j]).M() > (lepton0 + lepton1).M()){ // note that there could be negative invariant mass event
                DY_index = history_intersection.at(0); // set DY index when leptons selected
                if(verbose)
                    cout << "set lepton!" << endl;
                
                if(leptons[i]->Pt()>leptons[j]->Pt()){
                    lepton0 = *leptons[i];
                    lepton1 = *leptons[j];
                    
                    lepton0.SetIndexPIDStatus(leptons[i]->Index(), leptons[i]->PID(), leptons[i]->Status());
                    lepton1.SetIndexPIDStatus(leptons[j]->Index(), leptons[j]->PID(), leptons[j]->Status());
                }
                else{
                    lepton0 = *leptons[j];
                    lepton1 = *leptons[i];
                    lepton0.SetIndexPIDStatus(leptons[j]->Index(), leptons[j]->PID(), leptons[j]->Status());
                    lepton1.SetIndexPIDStatus(leptons[i]->Index(), leptons[i]->PID(), leptons[i]->Status());
                }
            }
        }// leptons(j) loop
    }// leptons(i) loop
    
    return DY_index;
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
    // [tunfold_hist]_[dipt-dimass]_[reco__fine_O-window_v1_UO]
    // [tunfold_hist]_[dipt-dimass]_[gen_dRp1__fine_O-window_v1_UO]
    // [tunfold_hist]_[dipt-dimass]_[gen_acceptance__fine_O-window_v1_UO]
    //
    // dipt_[reco__fine_O]_dimass_55to64
    
    // create 2d folded bins
    create_2d_folded_bin("dipt", "fine", false, true, "dimass", "window_v1", true, true);
    create_2d_folded_bin("dimass", "fine", true, true, "dipt", "window_v1", false, true);
    
    // create 2d unfolded bins
    create_2d_unfolded_bin("dipt", "coarse", false, true, "dimass", "window_v1", true, true);
    create_2d_unfolded_bin("dipt", "fine", false, true, "dimass", "window_v1", true, true);
    create_2d_unfolded_bin("dimass", "coarse", true, true, "dipt", "window_v1", false, true);
}
ISRAnalyzer::~ISRAnalyzer(){
    //DeleteCosThetaWeight();
}
