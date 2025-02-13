#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrectorWrapper.h"

#include <iostream>

void testFactorizedJetCorrectorWrapper() {
  
  FactorizedJetCorrectorWrapper *jec = new FactorizedJetCorrectorWrapper();
  jec->addJECset("Prompt24_V8M_DATA");

  jec->setRun(385900);
  jec->setJetEta(0.65);
  jec->setJetPhi(1.57);
  jec->setJetPt(130.);
  jec->setRho(25.);
  cout << "jec->getCorrection()=" << jec->getCorrection() << endl;

  jec->setRun(380900);
  jec->setJetEta(0.65);
  jec->setJetPhi(1.57);
  jec->setJetPt(130.);
  jec->setRho(25.);
  cout << "jec->getCorrection()=" << jec->getCorrection() << endl;

  cout << "Repeat jec->getCorrection()=" << jec->getCorrection() << endl;
}
