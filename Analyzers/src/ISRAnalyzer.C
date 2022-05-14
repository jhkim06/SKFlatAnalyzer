#include "ISRAnalyzer.h"

void ISRAnalyzer::initializeAnalyzer(){

    SMPAnalyzerCore::initializeAnalyzer(); // Z pt, Rochester, Z0
    
    if(HasFlag("nobjet")){
      vector<JetTagging::Parameters> jtps={JetTagging::Parameters(JetTagging::DeepCSV,JetTagging::Medium,JetTagging::mujets,JetTagging::mujets)};
      mcCorr->SetJetTaggingParameters(jtps);
    }
    
    IsNominalRun = !HasFlag("SYS") && !HasFlag("PDFSYS");
    
    if(fChain->GetListOfFiles()->GetEntries()){
      TString filename = fChain->GetListOfFiles()->At(0)->GetTitle();
      if(filename.Contains("SkimTree_")) IsSkimmed=true;
      else IsSkimmed = false;
    }else{
      cout<<"[ISRAnalyzer::initializeAnalyzer] no input file"<<endl;
      exit(EXIT_FAILURE);
    }
}

void ISRAnalyzer::executeEvent(){

    event=GetEvent();
    GetEventWeights(); // NEED TO STUDY,

    if(IsDYSample){

        vector<LHE> lhes=GetLHEs();
        LHE lhe_l0,lhe_l1;
        GetDYLHEParticles(lhes,lhe_l0,lhe_l1);

        const vector<Gen> gens=GetGens();
        Gen gen_parton0,gen_parton1,gen_l0,gen_l1,gen_l0_dressed,gen_l1_dressed,gen_l0_bare,gen_l1_bare;
        SMPAnalyzerCore::GetDYGenParticles(gens,gen_parton0,gen_parton1,gen_l0,gen_l1, PreFSR);
        SMPAnalyzerCore::GetDYGenParticles(gens,gen_parton0,gen_parton1,gen_l0_dressed,gen_l1_dressed, DressedDRp1);
        SMPAnalyzerCore::GetDYGenParticles(gens,gen_parton0,gen_parton1,gen_l0_bare,gen_l1_bare, PostFSR);

        Gen gen_isr_parton0, gen_isr_parton1, gen_isr_l0, gen_isr_l1, gen_isr_l0_bare, gen_isr_l1_bare;
        vector<const Gen*> added_photons;
        int DY_index = get_DY_gen_particles(gens, gen_isr_parton0, gen_isr_parton1, gen_isr_l0_bare, gen_isr_l1_bare, PostFSR);
        DY_index = get_DY_gen_particles(gens, gen_isr_parton0, gen_isr_parton1, gen_isr_l0, gen_isr_l1, PreFSR, added_photons);
        
        if(DY_index == -1) // DY leptons not properly selected
            return;
        
        /*
        cout << "hs dimass: " << (gen_l0_bare+gen_l1_bare).M() << " jh dimass: " << (gen_isr_l0_bare+gen_isr_l1_bare).M() << " DY index: " << DY_index << endl;
        cout << "hs index 0: " << gen_l0_bare.Index() << " hs index 1: " << gen_l1_bare.Index() << endl;
        cout << "jh index 0: " << gen_isr_l0_bare.Index() << " jh index 1: " << gen_isr_l1_bare.Index() << endl;
        cout << "hs dimass: " << (gen_l0+gen_l1).M() << " jh dimass: " << (gen_isr_l0+gen_isr_l1).M() << " DY index: " << DY_index << endl;
        cout << "jh index 0: " << gen_isr_l0.Index() << " jh index 1: " << gen_isr_l1.Index() << endl;
        */
        
        if(abs(gen_isr_l0_bare.PID()) == MUON){
            FillHist("Muon_DiLepton_Mass_preFSR", (gen_isr_l0 + gen_isr_l1).M(), lumiweight, 3000, 0, 3000);
            FillHist("Muon_DiLepton_Mass_postFSR", (gen_isr_l0_bare + gen_isr_l1_bare).M(), lumiweight, 3000, 0, 3000);
            FillHist("Muon_DiLepton_Mass_preFSR_HS", (gen_l0 + gen_l1).M(), lumiweight, 3000, 0, 3000);
            FillHist("Muon_DiLepton_Mass_postFSR_HS", (gen_l0_bare + gen_l1_bare).M(), lumiweight, 3000, 0, 3000);
            for(auto photon: added_photons){
                if(gen_l0_bare.DeltaR(*photon) < gen_l1_bare.DeltaR(*photon))
                    FillHist("Muon_DR_added_photon", gen_l0_bare.DeltaR(*photon), lumiweight, 100, 0, 1.);
                else
                    FillHist("Muon_DR_added_photon", gen_l1_bare.DeltaR(*photon), lumiweight, 100, 0, 1.);
            }
        }
        if(abs(gen_isr_l0_bare.PID()) == ELECTRON){
            FillHist("Electron_DiLepton_Mass_preFSR", (gen_isr_l0 + gen_isr_l1).M(), lumiweight, 3000, 0, 3000);
            FillHist("Electron_DiLepton_Mass_postFSR", (gen_isr_l0_bare + gen_isr_l1_bare).M(), lumiweight, 3000, 0, 3000);
            FillHist("Electron_DiLepton_Mass_preFSR_HS", (gen_l0 + gen_l1).M(), lumiweight, 3000, 0, 3000);
            FillHist("Electron_DiLepton_Mass_postFSR_HS", (gen_l0_bare + gen_l1_bare).M(), lumiweight, 3000, 0, 3000);
            for(auto photon: added_photons){
                if(gen_l0_bare.DeltaR(*photon) < gen_l1_bare.DeltaR(*photon))
                    FillHist("Electron_DR_added_photon", gen_l0_bare.DeltaR(*photon), lumiweight, 100, 0, 1.);
                else
                    FillHist("Electron_DR_added_photon", gen_l1_bare.DeltaR(*photon), lumiweight, 100, 0, 1.);
            }
        }

      /*
        if (abs((gen_l0_bare+gen_l1_bare).M()-(gen_isr_l0_bare+gen_isr_l1_bare).M()) > 1e-5){
            
            //get_DY_gen_particles(gens, gen_isr_parton0, gen_isr_parton1, gen_isr_l0, gen_isr_l1, PreFSR);
            get_DY_gen_particles(gens, gen_isr_parton0, gen_isr_parton1, gen_isr_l0, gen_isr_l1, PreFSR, added_photons);
            cout << "hs pre fsr mass: " << (gen_l0+gen_l1).M() << " jh pre fsr mass: " << (gen_isr_l0+gen_isr_l1).M() << endl;
            cout << "selected photons \n";
            for(auto i: added_photons){
                cout << i->Index() << " ";
            }
            cout << "\n";

            print_gen_particles(gens);
        }
        */

        // lumiweight
    }

    //
    if(!PassMETFilter()) return;
    
    TString prefix="";
    if(HasFlag("bjet")) prefix+="bjet/";
    
    if(IsNominalRun){
      FillCutflow(prefix+tauprefix+"cutflow","lumi",lumiweight);
      FillCutflow(prefix+tauprefix+"cutflow","PU",lumiweight*PUweight);
      FillCutflow(prefix+tauprefix+"cutflow","prefire",lumiweight*PUweight*prefireweight);
      FillCutflow(prefix+tauprefix+"cutflow","zpt",lumiweight*PUweight*prefireweight*zptweight);
      FillCutflow(prefix+tauprefix+"cutflow","z0",lumiweight*PUweight*prefireweight*zptweight*z0weight);
    }
    
    // NEED TO STUDY!
    int n_bjet=0;
    if(HasFlag("nobjet")){
      std::vector<Jet> jets=GetJets("tightLepVeto",30,2.7);
      std::sort(jets.begin(),jets.end(),PtComparing);

      JetTagging::Parameters jtp = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::mujets, JetTagging::mujets);
      for(const auto& jet:jets)
        if(mcCorr->IsBTagged_2a(jtp,jet))
            n_bjet++;
      
      if(HasFlag("nobjet")&&n_bjet) return;
      if(IsNominalRun) FillCutflow(prefix+tauprefix+"cutflow","BJetCut",lumiweight*PUweight*prefireweight*zptweight*z0weight);
    }
    
    if(DataYear==2016){
      vector<TString> muontrigger={
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
        "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      };
      vector<TString> emutrigger={
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
      };
      if(!HasFlag("emu") && event.PassTrigger(muontrigger))
        if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName(prefix+"mm2016");
      if(!HasFlag("emu") && event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))
        if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventWithChannelName(prefix+"ee2016");
      if(HasFlag("emu") && event.PassTrigger(emutrigger))
        if(!IsDATA||DataStream.Contains("MuonEG")) executeEventWithChannelName(prefix+"em2016");
    }else if(DataYear==2017){
      vector<TString> emutrigger={
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
      };
      if(!HasFlag("emu") && event.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"))
        if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName(prefix+"mm2017");
      if(!HasFlag("emu") && event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
        if(!IsDATA||DataStream.Contains("DoubleEG")) executeEventWithChannelName(prefix+"ee2017");
      if(HasFlag("emu") && event.PassTrigger(emutrigger))
        if(!IsDATA||DataStream.Contains("MuonEG")) executeEventWithChannelName(prefix+"em2017");
    }else if(DataYear==2018){
      vector<TString> emutrigger={
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
      };
      if(!HasFlag("emu") && event.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"))
        if(!IsDATA||DataStream.Contains("DoubleMuon")) executeEventWithChannelName(prefix+"mm2018");
      if(!HasFlag("emu") && event.PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"))
        if(!IsDATA||DataStream.Contains("EGamma")) executeEventWithChannelName(prefix+"ee2018");
      if(HasFlag("emu")  && event.PassTrigger(emutrigger))
        if(!IsDATA||DataStream.Contains("MuonEG")) executeEventWithChannelName(prefix+"em2018");
    }
}

void ISRAnalyzer::executeEventWithChannelName(TString channelname){
    
  map<TString, vector<Muon>> map_muons;
  map<TString, vector<Electron>> map_electrons;
  map<TString, ISRParameter> map_parameter; // Parameter defined in SMPAnalyzer
  
  if(channelname.Contains(TRegexp("mm20[0-9][0-9]"))){
      
    ISRParameter p("IDISO_SF_MediumID_trkIsoLoose_Q","",{"Mu17Leg1_MediumID_trkIsoLoose_Q","Mu8Leg2_MediumID_trkIsoLoose_Q"}, 20., 10.);
    
    map_muons[""]=MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0,3);
    map_parameter[""]=p.Clone(MakeLeptonPointerVector(map_muons[""]),
                  (IsNominalRun?NominalWeight:0)
                  +(HasFlag("SYS")&&!IsDATA?SystematicWeight:0)
                  +(HasFlag("PDFSYS")&&!IsDATA?PDFWeight:0)
                  );
    
    if(HasFlag("SYS")){
      map_muons["_scale_up"]=MuonMomentumCorrection(map_muons[""],+1);
      map_parameter["_scale_up"]=p.Clone(MakeLeptonPointerVector(map_muons["_scale_up"]));
                     
      map_muons["_scale_down"]=MuonMomentumCorrection(map_muons[""],-1);
      map_parameter["_scale_down"]=p.Clone(MakeLeptonPointerVector(map_muons["_scale_down"]));
      
      map_muons["_noroccor"]=MuonMomentumCorrection(map_muons[""],0,-1);
      map_parameter["_noroccor"]=p.Clone(MakeLeptonPointerVector(map_muons["_noroccor"]));
    }
  }else if(channelname.Contains(TRegexp("ee20[0-9][0-9]"))){
    ISRParameter p("ID_SF_MediumID_Q",{"Ele23Leg1_MediumID_Q","Ele12Leg2_MediumID_Q"},25.,15.);
    
    map_electrons["_noroccor"]=SMPGetElectrons("passMediumID",0.0,2.5);
    map_electrons[""]=ElectronEnergyCorrection(map_electrons["_noroccor"],0,0);
    map_parameter[""]=p.Clone(MakeLeptonPointerVector(map_electrons[""]),
                  (IsNominalRun?NominalWeight:0)
                  +(HasFlag("SYS")&&!IsDATA?SystematicWeight:0)
                  +(HasFlag("PDFSYS")&&!IsDATA?PDFWeight:0)
                  );

    if(HasFlag("SYS")){
      map_electrons["_scale_up"]=ScaleElectrons(map_electrons[""],1);
      std::sort(map_electrons["_scale_up"].begin(),map_electrons["_scale_up"].end(),PtComparing);
      map_parameter["_scale_up"]=p.Clone(MakeLeptonPointerVector(map_electrons["_scale_up"]));
      
      map_electrons["_scale_down"]=ScaleElectrons(map_electrons[""],-1);
      std::sort(map_electrons["_scale_down"].begin(),map_electrons["_scale_down"].end(),PtComparing);
      map_parameter["_scale_down"]=p.Clone(MakeLeptonPointerVector(map_electrons["_scale_down"]));
      
      map_electrons["_smear_up"]=SmearElectrons(map_electrons[""],1);
      std::sort(map_electrons["_smear_up"].begin(),map_electrons["_smear_up"].end(),PtComparing);
      map_parameter["_smear_up"]=p.Clone(MakeLeptonPointerVector(map_electrons["_smear_up"]));
      
      map_electrons["_smear_down"]=SmearElectrons(map_electrons[""],-1);
      std::sort(map_electrons["_smear_down"].begin(),map_electrons["_smear_down"].end(),PtComparing);
      map_parameter["_smear_down"]=p.Clone(MakeLeptonPointerVector(map_electrons["_smear_down"]));

      map_electrons["_eta2p5"]=ElectronEnergyCorrection(SMPGetElectrons("passMediumID",0.0,2.5),0,0);
      map_parameter["_eta2p5"]=p.Clone(MakeLeptonPointerVector(map_electrons["_eta2p5"]));

      map_parameter["_noroccor"]=p.Clone(MakeLeptonPointerVector(map_electrons["_noroccor"]));
      
      map_electrons["_noEcor"]=ElectronEnergyCorrection(map_electrons[""],-1,0);
      map_parameter["_noEcor"]=p.Clone(MakeLeptonPointerVector(map_electrons["_noEcor"]));

    }
  }else if(channelname.Contains(TRegexp("em20[0-9][0-9]"))){
    ISRParameter p;
    p.electronIDSF="ID_SF_MediumID_Q";
    p.muonIDSF="IDISO_SF_MediumID_trkIsoLoose_Q";
    p.triggerSF={"",""};
    p.lep0ptcut=25.;
    p.lep1ptcut=15.;

    map_electrons["_noroccor"]=SMPGetElectrons("passMediumID",0.0,2.4);
    map_electrons[""]=ElectronEnergyCorrection(map_electrons["_noroccor"],0,0);
    map_muons[""]=MuonMomentumCorrection(SMPGetMuons("POGMediumWithLooseTrkIso",0.0,2.4),0);

    std::vector<Lepton *> emu = MakeLeptonPointerVector(map_electrons[""]);
    std::vector<Lepton *> mu = MakeLeptonPointerVector(map_muons[""]);
    emu.insert(emu.end(),mu.begin(),mu.end());
    std::sort(emu.begin(),emu.end(),PtComparingPtr);

    map_parameter[""]=p.Clone(emu,
                  (IsNominalRun?NominalWeight:0)
                  +(HasFlag("SYS")&&!IsDATA?SystematicWeight:0)
                  +(HasFlag("PDFSYS")&&!IsDATA?PDFWeight:0)
                  );

  }else{
    cout<<"[ISRAnalyzer::executeEventWithPrefix] wrong channelname"<<endl;
    return;
  }
  
  ///////////////////////lepton selection///////////////////////
  for(const auto& [suffix,p]:map_parameter){
      TString prefix=tauprefix;
      double eventweight=lumiweight*PUweight*prefireweight*z0weight*zptweight;

      if(p.weightbit&NominalWeight) FillHist(channelname+"/"+prefix+"nlepton"+suffix,p.leps.size(),eventweight,10,0,10);
      
      if(p.leps.size()>=2){
          if(HasFlag("REGION_cf")){
              if(p.leps.at(0)->Charge()>0 && p.leps.at(1)->Charge()>0) prefix="pp_"+prefix;
              else if(p.leps.at(0)->Charge()<0 && p.leps.at(1)->Charge()<0) prefix="mm_"+prefix;
              else continue;
              
          }else{
              if(p.leps.at(0)->Charge() * p.leps.at(1)->Charge()>0) prefix="ss_"+prefix;
          }
          if(channelname.Contains(TRegexp("em20[0-9][0-9]"))){
              if(p.leps.at(0)->LeptonFlavour() == p.leps.at(1)->LeptonFlavour()) continue;
          }
          if(p.weightbit&NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"dilepton",eventweight);
          
          if(p.leps.at(0)->Pt() > p.lep0ptcut && p.leps.at(1)->Pt() > p.lep1ptcut){ // lepton eta cut applied in lepton selection,
              
              if(p.weightbit & NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"ptcut",eventweight);
              
              if(p.leps.at(0)->Charge() * p.leps.at(1)->Charge() < 0){
                  if(p.weightbit & NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"OS",eventweight);
                  
                 
                      /////////////////efficiency scale factors///////////////////
                      double IDSF=1.,IDSF_up=1.,IDSF_down=1.;
                      double ISOSF=1.,ISOSF_up=1.,ISOSF_down=1.;
                      double RECOSF=1.,RECOSF_up=1.,RECOSF_down=1.;
                      if(!IsDATA){
                          for(const auto& lep:p.leps){
                              TString LeptonIDSF_key="";
                              if(lep->LeptonFlavour()==Lepton::ELECTRON){
                                  LeptonIDSF_key=p.electronIDSF;

                                  double this_pt,this_eta;
                                  this_pt=((Electron*)lep)->UncorrPt();
                                  this_eta=((Electron*)lep)->scEta();
                                  
                                  double this_RECOSF=mcCorr->ElectronReco_SF(this_eta,this_pt,0);
                                  double this_RECOSF_up=mcCorr->ElectronReco_SF(this_eta,this_pt,1);
                                  double this_RECOSF_down=mcCorr->ElectronReco_SF(this_eta,this_pt,-1);
                                  RECOSF*=this_RECOSF; RECOSF_up*=this_RECOSF_up; RECOSF_down*=this_RECOSF_down;
                                  
                              }else if(lep->LeptonFlavour()==Lepton::MUON){
                                  LeptonIDSF_key=p.muonIDSF;
                                  double this_ISOSF=Lepton_SF(p.muonISOSF,lep,0);
                                  double this_ISOSF_up=Lepton_SF(p.muonISOSF,lep,1);
                                  double this_ISOSF_down=Lepton_SF(p.muonISOSF,lep,-1);
                                  ISOSF*=this_ISOSF; ISOSF_up*=this_ISOSF_up; ISOSF_down*=this_ISOSF_down;
                                  
                              }
                              double this_IDSF=Lepton_SF(LeptonIDSF_key,lep,0);
                              double this_IDSF_up=Lepton_SF(LeptonIDSF_key,lep,1);
                              double this_IDSF_down=Lepton_SF(LeptonIDSF_key,lep,-1);
                              IDSF*=this_IDSF; IDSF_up*=this_IDSF_up; IDSF_down*=this_IDSF_down;
                              
                          }
                  
                      }
      
                      double triggerSF=1.,triggerSF_up=1.,triggerSF_down=1.;
                      if(!IsDATA){
                          if(p.triggerSF.size()==1){
                              triggerSF*=LeptonTrigger_SF(p.triggerSF[0],p.leps,0);
                              triggerSF_up*=LeptonTrigger_SF(p.triggerSF[0],p.leps,1);
                              triggerSF_down*=LeptonTrigger_SF(p.triggerSF[0],p.leps,-1);
                              
                          }else if(p.triggerSF.size()==2){
                              triggerSF*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,0);
                              triggerSF_up*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,1);
                              triggerSF_down*=DileptonTrigger_SF(p.triggerSF[0],p.triggerSF[1],p.leps,-1);
                              
                          }
                          
                      }
                      if(p.weightbit&NominalWeight){
                          FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"RECO",eventweight*RECOSF);
                          FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"ID",eventweight*RECOSF*IDSF);
                          FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"ISO",eventweight*RECOSF*IDSF*ISOSF);
                          FillCutflow(channelname+"/"+prefix+"cutflow"+suffix,"trigger",eventweight*RECOSF*IDSF*ISOSF*triggerSF);
                      }
    
                      ///////////////////////map_weight//////////////////
                      map<TString,double> map_weight;
                      map<TString,double> map_gen_weight;
                      map<TString,double> map_reco_weight;
                      if(p.weightbit&NominalWeight){
                          map_weight[""]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_gen_weight[""]=lumiweight*PUweight*zptweight;
                          map_reco_weight[""]=prefireweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          
                      }
                      if(p.weightbit&SystematicWeight){
                          // only for the nominal parameter selection
                          map_weight["_noPUweight"]=lumiweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_weight["_PUweight_up"]=lumiweight*PUweight_up*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_weight["_PUweight_down"]=lumiweight*PUweight_down*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_gen_weight["_PUweight_up"]=lumiweight*PUweight_up*zptweight;
                          map_reco_weight["_PUweight_up"]=prefireweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_gen_weight["_PUweight_down"]=lumiweight*PUweight_down*zptweight;
                          map_reco_weight["_PUweight_down"]=prefireweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
              
                          map_weight["_noprefireweight"]=lumiweight*PUweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_weight["_prefireweight_up"]=lumiweight*PUweight*prefireweight_up*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_weight["_prefireweight_down"]=lumiweight*PUweight*prefireweight_down*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_gen_weight["_prefireweight_up"]=lumiweight*PUweight*zptweight;
                          map_reco_weight["_prefireweight_up"]=prefireweight_up*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_gen_weight["_prefireweight_down"]=lumiweight*PUweight*zptweight;
                          map_reco_weight["_prefireweight_down"]=prefireweight_down*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          
                          map_weight["_nozptweight"]=lumiweight*PUweight*prefireweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_weight["_noz0weight"]=lumiweight*PUweight*prefireweight*zptweight*RECOSF*IDSF*ISOSF*triggerSF;
                          map_weight["_noefficiencySF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight;
              
                          map_weight["_noRECOSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*IDSF*ISOSF*triggerSF;
                          map_weight["_RECOSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF_up*IDSF*ISOSF*triggerSF;
                          map_weight["_RECOSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF_down*IDSF*ISOSF*triggerSF;
                          map_gen_weight["_RECOSF_up"]=lumiweight*PUweight*zptweight;
                          map_reco_weight["_RECOSF_up"]=prefireweight*z0weight*RECOSF_up*IDSF*ISOSF*triggerSF;
                          map_gen_weight["_RECOSF_down"]=lumiweight*PUweight*zptweight;
                          map_reco_weight["_RECOSF_down"]=prefireweight*z0weight*RECOSF_down*IDSF*ISOSF*triggerSF;
              
                          map_weight["_noIDSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*ISOSF*triggerSF;
                          map_weight["_IDSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF_up*ISOSF*triggerSF;
                          map_weight["_IDSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF_down*ISOSF*triggerSF;
                          map_gen_weight["_IDSF_up"]=lumiweight*PUweight*zptweight;
                          map_reco_weight["_IDSF_up"]=prefireweight*z0weight*RECOSF*IDSF_up*ISOSF*triggerSF;
                          map_gen_weight["_IDSF_down"]=lumiweight*PUweight*zptweight;
                          map_reco_weight["_IDSF_down"]=prefireweight*z0weight*RECOSF*IDSF_down*ISOSF*triggerSF;
              
                          map_weight["_noISOSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*triggerSF;
                          map_weight["_ISOSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF_up*triggerSF;
                          map_weight["_ISOSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF_down*triggerSF;
                          map_gen_weight["_ISOSF_up"]=lumiweight*PUweight*zptweight;
                          map_reco_weight["_ISOSF_up"]=prefireweight*z0weight*RECOSF*IDSF*ISOSF_up*triggerSF;
                          map_gen_weight["_ISOSF_down"]=lumiweight*PUweight*zptweight;
                          map_reco_weight["_ISOSF_down"]=prefireweight*z0weight*RECOSF*IDSF*ISOSF_down*triggerSF;
              
                          map_weight["_notriggerSF"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF;
                          map_weight["_triggerSF_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF_up;
                          map_weight["_triggerSF_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF_down;
                          map_gen_weight["_triggerSF_up"]=lumiweight*PUweight*zptweight;
                          map_reco_weight["_triggerSF_up"]=prefireweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF_up;
                          map_gen_weight["_triggerSF_down"]=lumiweight*PUweight*zptweight;
                          map_reco_weight["_triggerSF_down"]=prefireweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF_down;
                          
                      }
                      if(p.weightbit&PDFWeight){
                          for(unsigned int i=0;i<PDFWeights_Scale->size();i++){
                              map_weight[Form("_scalevariation%d",i)]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_Scale->at(i);
                              
                          }
                          for(unsigned int i=0;i<PDFWeights_Error->size();i++){
                              map_weight[Form("_pdf%d",i)]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_Error->at(i);
                              
                          }
                          if(PDFWeights_AlphaS->size()==2){
                              map_weight["_alphaS_up"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_AlphaS->at(0);
                              map_weight["_alphaS_down"]=lumiweight*PUweight*prefireweight*zptweight*z0weight*RECOSF*IDSF*ISOSF*triggerSF*PDFWeights_AlphaS->at(1);
                          }
                          
                      }
                  
                      if ((*p.leps.at(0)+*p.leps.at(1)).Pt()<p.dilep_pt_cut){ // dilepton pt cut, (100 GeV default)
                          if(p.weightbit & NominalWeight) FillCutflow(channelname+"/"+prefix+"cutflow"+suffix, "dipt cut", eventweight);
                          
                          fill_unfold_hists(channelname, prefix, suffix, (Particle*)p.leps[0], (Particle*)p.leps[1], map_weight, *tunfold_parameter, TUnfold_Bin::smeared_bin);
                         
                          if (IsDYSample && prefix==""){ // prefix "" means non tautau event
                              // response matrix
                            
                              const vector<Gen> gens=GetGens();
                              Gen gen_isr_parton0, gen_isr_parton1, gen_isr_l0, gen_isr_l1, gen_isr_l0_bare, gen_isr_l1_bare;
                              vector<const Gen*> added_photons;
                              int DY_index = get_DY_gen_particles(gens, gen_isr_parton0, gen_isr_parton1, gen_isr_l0, gen_isr_l1, PreFSR, added_photons);
                              fill_unfold_hists(channelname, prefix, suffix, (Particle*)&gen_isr_l0, (Particle*)&gen_isr_l1, map_weight, *tunfold_parameter, TUnfold_Bin::truth_bin);
                              fill_unfold_response_matrixs(channelname, prefix, suffix, (Particle*)p.leps[0], (Particle*)p.leps[1], (Particle*)&gen_isr_l0, (Particle*)&gen_isr_l1,
                                                map_reco_weight, map_gen_weight, *tunfold_parameter);
                          }
                          
                      
                      // dilepton pt cut
                      /*
                       variable bin, pt cut not needed
                       */

                      ///////////////////////fill hists///////////////////////
                      if(HasFlag("TOY")){
                          //FillHistsToy(channelname,prefix,suffix,(Particle*)p.leps[0],(Particle*)p.leps[1],map_weight);
                      }
                      else{
                          fill_ISR_hists(channelname, prefix, suffix, (Particle*)p.leps[0], (Particle*)p.leps[1], map_weight);
                          if(IsDYSample&&prefix==""&&IsNominalRun){
                              vector<Gen> gens=GetGens();
                              Gen truth_l0=GetGenMatchedLepton(*p.leps[0],gens);
                              Gen truth_l1=GetGenMatchedLepton(*p.leps[1],gens);
                              if(!truth_l0.IsEmpty()&&!truth_l1.IsEmpty()){
                                  //FillHists(channelname,"truth_",suffix,(Particle*)&truth_l0,(Particle*)&truth_l1,map_weight);
                              }
                              //else cout<<"no matching"<<endl;
                          }
                      }
                      // effect of over/underflow bin
                      
                  }// dilepton pt cut
              } // OS
          } // lepton pt cut
      } // two lepton
  } // map_parameters
}

void ISRAnalyzer::fill_ISR_hists(TString channelname, TString pre, TString suf, Particle* l0, Particle* l1, map<TString,double> map_weight){
    
    // as first try, get only essential histogram.
    
    TLorentzVector dilepton = (*l0) + (*l1);
    double dimass = dilepton.M();
    double dipt = dilepton.Pt();
   
    // mass dependent histograms
    for(int i = 0; i < nmass_window; i++){
        double low_mass_edge = mass_window[i];
        double high_mass_edge = mass_window[i+1];
        
        double* mass_bin_pointer = (double*)mass_bin_fine_muon;
        int n_mass_bin = n_mass_bin_fine_muon-1;
        
        if(i == 0 && channelname.Contains(TRegexp("ee20[0-9][0-9]"))){
            low_mass_edge = 50;
            mass_bin_pointer = (double*)mass_bin_fine_muon;
            n_mass_bin = n_mass_bin_fine_electron-1;
        }
        
        if(dimass > low_mass_edge && dimass < high_mass_edge){
            
            string m = "m";
            string to = "to";
            string mass_window_postfix = m + Form("%d", (int)low_mass_edge) + to + Form("%d", (int)high_mass_edge);
            
            FillHist(channelname+"/"+pre+"dilep_pt_m"+ mass_window_postfix +suf, dipt, map_weight, sizeof(pt_bin)/sizeof(double)-1, (double*)pt_bin);
            FillHist(channelname+"/"+pre+"dilep_mass_m"+ mass_window_postfix +suf, dimass, map_weight, n_mass_bin, mass_bin_pointer);
        }
    }// loop mass window
}

void ISRAnalyzer::executeEventFromParameter(AnalyzerParameter param){

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

int ISRAnalyzer::get_DY_gen_particles(const vector<Gen>& gens, Gen& parton0, Gen& parton1, Gen& lepton0, Gen& lepton1, int mode){


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
    
    Gen lepton0_temp = lepton0;
    Gen lepton1_temp = lepton1;
    
    // loop leptons
    int nlepton = leptons.size();
    for(int i = 0;i < nlepton;i++){
        if(leptons[i]->Index() == lepton0.Index() || leptons[i]->Index() == lepton1.Index()) continue;
        for(int j = i + 1;j < nlepton;j++){
            if(leptons[j]->Index() == lepton0.Index()||leptons[j]->Index() == lepton1.Index()) continue;
            if(!(leptons[i]->PID()+leptons[j]->PID() == 0)) continue;
            vector<int> history_i = TrackGenSelfHistory(*leptons[i], gens);
            vector<int> history_j = TrackGenSelfHistory(*leptons[j], gens);
            if(history_i.at(1) == history_j.at(1)) photons.push_back(&gens[history_i.at(1)]);

        }
    }
    
    for(const auto photon : photons){

        if(mode == DressedMode::MotherMatch || mode == DressedMode::MotherDRMatch){
            auto it = find(DY_history.begin(), DY_history.end(), photon->MotherIndex());
            // photon's mother not exist in DY history
            if(it == DY_history.end()){
                continue;
            }
            
            // photon's mother is not lepton
            if(abs(gens.at(*it).PID()) != abs(lepton0_temp.PID())){
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
// TODO return DY vertex index

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

            // check if lepton vertex included
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
   
    tunfold_parameter = new TUnfoldParameter{sizeof(pt_bin_fine)/sizeof(double)-1, pt_bin_fine, sizeof(pt_bin_coarse)/sizeof(double)-1, pt_bin_coarse,
        sizeof(mass_window)/sizeof(double)-1, mass_window, sizeof(mass_window)/sizeof(double)-1, mass_window, false, true, true, true, "nominal_pt", "dipt", "dimass"};
    /*
    create_tunfold_hist(sizeof(pt_bin_fine)/sizeof(double)-1, pt_bin_fine,
                        sizeof(pt_bin_coarse)/sizeof(double)-1, pt_bin_coarse,
                        sizeof(mass_window)/sizeof(double)-1, mass_window,
                        false, true, true, true, "dipt", "dimass");
     */

}

ISRAnalyzer::~ISRAnalyzer(){

}


