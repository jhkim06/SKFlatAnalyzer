#include "ISRAnalyzer.h"

void ISRAnalyzer::initializeAnalyzer(){

}

void ISRAnalyzer::executeEvent(){


  AnalyzerParameter param;

  executeEventFromParameter(param);

}

void ISRAnalyzer::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

ISRAnalyzer::ISRAnalyzer(){

}

ISRAnalyzer::~ISRAnalyzer(){

}


