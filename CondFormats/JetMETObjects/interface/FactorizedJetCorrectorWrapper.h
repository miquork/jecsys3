// Wrapper for FactorizeJetCorrection to implement run-number based IOVs

#ifndef FACTORIZED_JET_CORRECTOR_WRAPPER_H
#define FACTORIZED_JET_CORRECTOR_WRAPPER_H

#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <map>

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

class FactorizedJetCorrectorWrapper
{
 public:
  FactorizedJetCorrectorWrapper() :
    mJetEta(0), mJetPhi(0), mJetPt(0), mJetA(0), mRho(0),
    mIsJetEtaset(false), mIsJetPhiset(false), mIsJetPtset(false),
    mIsJetAset(false), mIsRhoset(false), mIsRunset(false), mIsJECset(false),
    mCurrentRunMin(-1), mCurrentRunMax(-1), mJEC(0) {};
  void addJEC(FactorizedJetCorrector *fJEC, int fRunMin, int fRunMax) {
    mJECmap[std::make_pair(fRunMin,fRunMax)] = fJEC;
  }
  virtual void addJECset(std::string id, int verbosity = 2);
  void setJetEta      (float fEta) { mJetEta = fEta; mIsJetEtaset = true; }
  void setJetPhi      (float fPhi) { mJetPhi = fPhi; mIsJetPhiset = true; } 
  void setJetPt       (float fPt)  {   mJetPt = fPt;  mIsJetPtset = true; } 
  void setJetA        (float fA)   {     mJetA = fA;   mIsJetAset = true; } 
  void setRho         (float fRho) {    mRho = fRho;    mIsRhoset = true; }
  void setRun         (int fRun);
  float getCorrection();
  std::vector<float> getSubCorrections();
  
 private:
  //---- Member Data ---------
  float mJetEta;
  float mJetPhi;
  float mJetPt;
  float mJetA;
  float mRho;
  int   mRun;
  bool  mIsJetEtaset;
  bool  mIsJetPhiset;
  bool  mIsJetPtset;
  bool  mIsJetAset;
  bool  mIsRhoset;
  bool  mIsRunset;
  bool  mIsJECset;
  // ---- Cache ----
  int mCurrentRunMin;
  int mCurrentRunMax;
  FactorizedJetCorrector *mJEC;
  std::map<std::pair<int, int>, FactorizedJetCorrector*> mJECmap;
  typedef std::map<std::pair<int, int>, FactorizedJetCorrector*>::const_iterator IT;

  bool setJEC();
};

void FactorizedJetCorrectorWrapper::setRun(int fRun) {

  // Find correct FactorizedJetCorrector, if not cached already
  if (fRun<mCurrentRunMin || fRun>mCurrentRunMax) {
    mIsJECset = false;
    for (IT it = mJECmap.begin(); it != mJECmap.end(); ++it) {
      if (fRun>=it->first.first && fRun<=it->first.second) {
	mCurrentRunMin = it->first.first;
	mCurrentRunMax = it->first.second;
	mJEC = it->second;
	if (mJEC!=0) mIsJECset = true;
	continue;
      }
    }
  }
  mRun = fRun;
  mIsRunset = true;
} // setRun

// Pass variables to FactorizedJetCorrector and mark them as consumed
// Return true if successful
bool FactorizedJetCorrectorWrapper::setJEC() {

  if (mIsRunset && mIsJECset) {
    if (mIsRhoset)    mJEC->setRho(mRho);
    if (mIsJetAset)   mJEC->setJetA(mJetA);
    if (mIsJetPtset)  mJEC->setJetPt(mJetPt);
    if (mIsJetEtaset) mJEC->setJetEta(mJetEta);
    if (mIsJetPhiset) mJEC->setJetPhi(mJetPhi);
    mIsRhoset = mIsJetAset = mIsJetPtset = mIsJetEtaset = mIsJetPhiset = false;
    mIsRunset = mIsJECset = false;
    return true;
  }
  else
    return false;
} // setJEC

float FactorizedJetCorrectorWrapper::getCorrection() {

  if (setJEC()) return mJEC->getCorrection();
  else {
    std::cerr << "ERROR: inputs not properly set. Exiting" << endl;
    exit(0);
  }
} // getCorrection

std::vector<float> FactorizedJetCorrectorWrapper::getSubCorrections() {

  if (setJEC()) return mJEC->getSubCorrections();
  else {
    std::cerr << "ERROR: inputs not properly set. Exiting" << endl;
    exit(0);
  }
} // getSubCorrection

#endif
