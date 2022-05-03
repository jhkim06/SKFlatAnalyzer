#ifndef ISRAnalyzer_h
#define ISRAnalyzer_h

#include "AnalyzerCore.h"

class ISRAnalyzer : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  ISRAnalyzer();
  ~ISRAnalyzer();

};



#endif

