#include "TelEventTTreeReader.hpp"
#include <TFile.h>
#include <TTree.h>

#include <iostream>

void altel::TelEventTTreeReader::setTTree(TTree *pTTree){
  if(!pTTree){
    std::fprintf(stderr, "TTree is not yet set\n");
    throw;
  }
  m_pTTree = pTTree;
  TTree &tree = *m_pTTree;
  m_numEvents = tree.GetEntries();
  tree.ResetBranchAddresses();

  tree.SetBranchAddress("RunN", &rRunN);
  tree.SetBranchAddress("EventN", &rEventN);
  tree.SetBranchAddress("ConfigN", &rConfigN);
  tree.SetBranchAddress("Clock", (ULong64_t*)&rClock);
  tree.SetBranchAddress("NumTraj_PerEvent", &rNumTraj_PerEvent);

  tree.SetBranchAddress("RawMeasVec_DetN", &pRawMeasVec_DetN);
  tree.SetBranchAddress("RawMeasVec_U", &pRawMeasVec_U);
  tree.SetBranchAddress("RawMeasVec_V", &pRawMeasVec_V);
  tree.SetBranchAddress("RawMeasVec_Clk", &pRawMeasVec_Clk);

  tree.SetBranchAddress("HitMeasVec_DetN", &pHitMeasVec_DetN);
  tree.SetBranchAddress("HitMeasVec_U", &pHitMeasVec_U);
  tree.SetBranchAddress("HitMeasVec_V", &pHitMeasVec_V);
  tree.SetBranchAddress("HitMeasVec_NumRawMeas_PerHitMeas", &pHitMeasVec_NumRawMeas_PerHitMeas);
  tree.SetBranchAddress("HitMeasVec_Index_To_RawMeas", &pHitMeasVec_Index_To_RawMeas);

  tree.SetBranchAddress("HitFitVec_DetN", &pHitFitVec_DetN);
  tree.SetBranchAddress("HitFitVec_U", &pHitFitVec_U);
  tree.SetBranchAddress("HitFitVec_V", &pHitFitVec_V);
  tree.SetBranchAddress("HitFitVec_X", &pHitFitVec_X);
  tree.SetBranchAddress("HitFitVec_Y", &pHitFitVec_Y);
  tree.SetBranchAddress("HitFitVec_Z", &pHitFitVec_Z);
  tree.SetBranchAddress("HitFitVec_DirX", &pHitFitVec_DirX);
  tree.SetBranchAddress("HitFitVec_DirY", &pHitFitVec_DirY);
  tree.SetBranchAddress("HitFitVec_DirZ", &pHitFitVec_DirZ);
  tree.SetBranchAddress("HitFitVec_Index_To_Origin_HitMeas", &pHitFitVec_Index_To_Origin_HitMeas);
  tree.SetBranchAddress("HitFitVec_Index_To_Matched_HitMeas", &pHitFitVec_Index_To_Matched_HitMeas);

  tree.SetBranchAddress("TrajVec_NumHitFit_PerTraj", &pTrajVec_NumHitFit_PerTraj);
  tree.SetBranchAddress("TrajVec_NumHitMeas_Origin_PerTraj", &pTrajVec_NumHitMeas_Origin_PerTraj);
  tree.SetBranchAddress("TrajVec_NumHitMeas_Matched_PerTraj", &pTrajVec_NumHitMeas_Matched_PerTraj);

  tree.SetBranchAddress("TrajVec_Index_To_HitFit", &pTrajVec_Index_To_HitFit);

  // ana
  tree.SetBranchAddress("AnaVec_Matched_DetN", &pAnaVec_Matched_DetN);
  tree.SetBranchAddress("AnaVec_Matched_ResdU", &pAnaVec_Matched_ResdU);
  tree.SetBranchAddress("AnaVec_Matched_ResdV", &pAnaVec_Matched_ResdV);
}

size_t altel::TelEventTTreeReader::numEvents() const{
  return m_numEvents;
}

std::shared_ptr<altel::TelEvent> altel::TelEventTTreeReader::createTelEvent(size_t n){
  if(n>=m_numEvents){
    std::fprintf(stderr, "no more entries in ttree");
    return nullptr;
  }

  m_pTTree->GetEntry(n);
  std::shared_ptr<altel::TelEvent> telEvent(new altel::TelEvent(rRunN, rEventN, rConfigN, rClock));

  size_t numMeasHit = rHitMeasVec_DetN.size();
  assert(rHitMeasVec_U.size() == numMeasHit && rHitMeasVec_V.size() == numMeasHit);
  std::vector<std::shared_ptr<altel::TelMeasHit>> measHits;
  measHits.reserve(rHitMeasVec_DetN.size());

  auto it_measHitVec_detN = rHitMeasVec_DetN.begin();
  auto it_measHitVec_U = rHitMeasVec_U.begin();
  auto it_measHitVec_V = rHitMeasVec_V.begin();
  auto it_measHitVec_detN_end = rHitMeasVec_DetN.end();
  while(it_measHitVec_detN !=it_measHitVec_detN_end){
    measHits.emplace_back(new altel::TelMeasHit(*it_measHitVec_detN,
                                                *it_measHitVec_U,
                                                *it_measHitVec_V,
                                                {}));
    it_measHitVec_detN++;
    it_measHitVec_U++;
    it_measHitVec_V++;
  }

  std::vector<std::shared_ptr<altel::TelTrajHit>> trajHits;
  size_t numFitHit = rHitFitVec_DetN.size();
  trajHits.reserve(numFitHit);

  auto it_fitHitVec_detN = rHitFitVec_DetN.begin();
  auto it_fitHitVec_U = rHitFitVec_U.begin();
  auto it_fitHitVec_V = rHitFitVec_V.begin();
  auto it_fitHitVec_X = rHitFitVec_X.begin();
  auto it_fitHitVec_Y = rHitFitVec_Y.begin();
  auto it_fitHitVec_Z = rHitFitVec_Z.begin();
  auto it_fitHitVec_DX = rHitFitVec_DirX.begin();
  auto it_fitHitVec_DY = rHitFitVec_DirY.begin();
  auto it_fitHitVec_DZ = rHitFitVec_DirZ.begin();
  auto it_fitHitVec_OriginMeasHit_index = rHitFitVec_Index_To_Origin_HitMeas.begin();
  auto it_fitHitVec_MatchedMeasHit_index = rHitFitVec_Index_To_Matched_HitMeas.begin();

  auto it_fitHitVec_detN_end = rHitFitVec_DetN.end();

  while(it_fitHitVec_detN !=it_fitHitVec_detN_end){
    uint16_t originMeasHitIndex = *it_fitHitVec_OriginMeasHit_index;
    std::shared_ptr<altel::TelMeasHit> originMeasHit;
    if(originMeasHitIndex!=int16_t(-1)){
      assert( originMeasHitIndex<measHits.size() );
      originMeasHit = measHits[originMeasHitIndex];
    }

    auto fitHit = std::make_shared<altel::TelFitHit>(*it_fitHitVec_detN,
                                                     *it_fitHitVec_U,
                                                     *it_fitHitVec_V,
                                                     *it_fitHitVec_X,
                                                     *it_fitHitVec_Y,
                                                     *it_fitHitVec_Z,
                                                     *it_fitHitVec_DX,
                                                     *it_fitHitVec_DY,
                                                     *it_fitHitVec_DZ,
                                                     originMeasHit);


    uint16_t matchedMeasHitIndex = *it_fitHitVec_MatchedMeasHit_index;
    std::shared_ptr<altel::TelMeasHit> matchedMeasHit;
    if(matchedMeasHitIndex!= int16_t(-1)){
      assert( matchedMeasHitIndex<measHits.size() );
      matchedMeasHit = measHits[matchedMeasHitIndex];
    }

    trajHits.emplace_back(new altel::TelTrajHit(*it_fitHitVec_detN,
                                                fitHit,
                                                matchedMeasHit));

    it_fitHitVec_detN ++;
    it_fitHitVec_U ++;
    it_fitHitVec_V ++;
    it_fitHitVec_X ++;
    it_fitHitVec_Y ++;
    it_fitHitVec_Z ++;
    it_fitHitVec_DX ++;
    it_fitHitVec_DY ++;
    it_fitHitVec_DZ ++;
    it_fitHitVec_OriginMeasHit_index++;
    it_fitHitVec_MatchedMeasHit_index++;
  }

  std::vector<std::shared_ptr<altel::TelTrajectory>> trajs;

  auto it_numFitHit_PerTraj = rTrajVec_NumHitFit_PerTraj.begin();
  auto it_numFitHit_PerTraj_end = rTrajVec_NumHitFit_PerTraj.end();
  auto it_fitHit_index = rTrajVec_Index_To_HitFit.begin();
  auto it_fitHit_index_end = rTrajVec_Index_To_HitFit.end();

  while(it_numFitHit_PerTraj != it_numFitHit_PerTraj_end){
    auto traj = std::make_shared<altel::TelTrajectory>();
    size_t numFitHits_got = 0;
    while(it_fitHit_index!=it_fitHit_index_end && numFitHits_got < *it_numFitHit_PerTraj){
      int16_t fitHitIndex = *it_fitHit_index;
      if(fitHitIndex!=int16_t(-1)){
        auto trajHit = trajHits[fitHitIndex];
        traj->trajectoryHits().push_back(trajHit);
        numFitHits_got++;
      }
      it_fitHit_index++;
    }
    trajs.push_back(traj);
    it_numFitHit_PerTraj++;
  }

  telEvent->measHits() = std::move(measHits);
  telEvent->trajectories() = std::move(trajs);

  return telEvent;
}
