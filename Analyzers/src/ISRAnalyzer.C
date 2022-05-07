#include "ISRAnalyzer.h"

void ISRAnalyzer::initializeAnalyzer(){

    SMPAnalyzerCore::initializeAnalyzer();
}

void ISRAnalyzer::executeEvent(){

    event=GetEvent();
    GetEventWeights();

    if(IsDYSample){

        /*
         dilepton pt/mass distribution
         at lhe
         at pre-fsr
         at post-fsr


         isPromptFinalState: prompt and final state(status1)
         require the mother is not proton

         post-fsr
         loop over all lepton pairs
         find a pair which is not from lepton vertex,

         pre-fsr
         photon's mother is in the track history


         SMPAnalyzerCore::GetDYGenParticles()
         post-fsr
         - select prompt and status1 leptons
         - find same flavour, opposite sign pairs with largest mass
         pre-fsr
         - if PID of photon's mother is the same with that of post-fsr lepton add it.
         - if the remaining pair's mother is the same think is as photon
         -
         */

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

        
        if (abs((gen_l0_bare+gen_l1_bare).M()-(gen_isr_l0_bare+gen_isr_l1_bare).M()) > 1e-5){
            cout << "hs dimass: " << (gen_l0_bare+gen_l1_bare).M() << " jh dimass: " << (gen_isr_l0_bare+gen_isr_l1_bare).M() << " DY index: " << DY_index << endl;
            cout << "hs index 0: " << gen_l0_bare.Index() << " hs index 1: " << gen_l1_bare.Index() << endl;
            cout << "jh index 0: " << gen_isr_l0_bare.Index() << " jh index 1: " << gen_isr_l1_bare.Index() << endl;
            
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
       

        // lumiweight
        // function of mass bins
        // dilep mass, dilepton pt, lepton pt, lepton eta,

        // Compare with my Gen particle selection

    }

    //
    if(!PassMETFilter()) return;

    //AnalyzerParameter param;
    //executeEventFromParameter(param);
}

void ISRAnalyzer::executeEventFromParameter(AnalyzerParameter param){
    //if(!PassMETFilter()) return;
    //Event ev = GetEvent();

}

void ISRAnalyzer::FillHists(const TString channelname, const TString pre, const TString suf,
                            const Particle* l0, const Particle* l1, const double weight) const {


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

  
    DY_index = get_DY_bare_lepton_pair(gens, leptons, lepton0, lepton1);


    return DY_index;
}

//
int ISRAnalyzer::get_DY_dressed_lepton_pair(const vector<Gen>& gens, const vector<const Gen*>& leptons, vector<const Gen*>& photons, Gen& lepton0, Gen& lepton1,
                                            const DressedMode mode, vector<const Gen*>& added_photons, const double dR){

    
    const int DY_index = get_DY_bare_lepton_pair(gens, leptons, lepton0, lepton1);
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
    
    //TODO save features of the added FSR gamma

    return DY_index;
}


int ISRAnalyzer::get_DY_bare_lepton_pair(const vector<Gen>& gens, const vector<const Gen*>& leptons, Gen& lepton0, Gen& lepton1){
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
            DY_index = history_intersection.at(0);

            // set DY lepton pair
            if((*leptons[i] + *leptons[j]).M() > (lepton0 + lepton1).M()){
                if(leptons[i]->Pt()>leptons[j]->Pt()){
                    lepton0 = *leptons[i];
                    lepton1 = *leptons[j];
                }
                else{
                    lepton0 = *leptons[j];
                    lepton1 = *leptons[i];
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

}

ISRAnalyzer::~ISRAnalyzer(){

}


