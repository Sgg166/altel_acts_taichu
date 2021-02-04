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
  tree.SetBranchAddress("EveN", &rEventN);
  tree.SetBranchAddress("DetN", &rConfigN);
  tree.SetBranchAddress("ClkN", (ULong64_t*)&rClock);
  tree.SetBranchAddress("NumTrajs_PerEvent", &rNumTraj_PerEvent);
  tree.SetBranchAddress("NumMeasHits_PerEvent", &rNumMeasHits_PerEvent);

  tree.SetBranchAddress("MeasRawVec_DetN", &pRawMeasVec_DetN);
  tree.SetBranchAddress("MeasRawVec_U", &pRawMeasVec_U);
  tree.SetBranchAddress("MeasRawVec_V", &pRawMeasVec_V);
  tree.SetBranchAddress("MeasRawVec_Clk", &pRawMeasVec_Clk);

  tree.SetBranchAddress("MeasHitVec_DetN", &pHitMeasVec_DetN);
  tree.SetBranchAddress("MeasHitVec_U", &pHitMeasVec_U);
  tree.SetBranchAddress("MeasHitVec_V", &pHitMeasVec_V);
  tree.SetBranchAddress("MeasHitVec_NumMeasRaws_PerMeasHit", &pHitMeasVec_NumRawMeas_PerHitMeas);
  tree.SetBranchAddress("MeasHitVec_MeasRaw_Index", &pHitMeasVec_Index_To_RawMeas);

  tree.SetBranchAddress("TrajHitVec_DetN", &pHitFitVec_DetN);
  tree.SetBranchAddress("TrajHitVec_U", &pHitFitVec_U);
  tree.SetBranchAddress("TrajHitVec_V", &pHitFitVec_V);
  tree.SetBranchAddress("TrajHitVec_X", &pHitFitVec_X);
  tree.SetBranchAddress("TrajHitVec_Y", &pHitFitVec_Y);
  tree.SetBranchAddress("TrajHitVec_Z", &pHitFitVec_Z);
  tree.SetBranchAddress("TrajHitVec_DirX", &pHitFitVec_DirX);
  tree.SetBranchAddress("TrajHitVec_DirY", &pHitFitVec_DirY);
  tree.SetBranchAddress("TrajHitVec_DirZ", &pHitFitVec_DirZ);
  tree.SetBranchAddress("TrajHitVec_OriginMeasHit_Index", &pHitFitVec_Index_To_Origin_HitMeas);
  tree.SetBranchAddress("TrajHitVec_MatchedMeasHit_Index", &pHitFitVec_Index_To_Matched_HitMeas);

  tree.SetBranchAddress("TrajVec_NumTrajHits_PerTraj", &pTrajVec_NumHitFit_PerTraj);
  tree.SetBranchAddress("TrajVec_NumOriginMeasHits_PerTraj", &pTrajVec_NumHitMeas_Origin_PerTraj);
  tree.SetBranchAddress("TrajVec_NumMatchedMeasHits_PerTraj", &pTrajVec_NumHitMeas_Matched_PerTraj);

  tree.SetBranchAddress("TrajVec_TrajHit_Index", &pTrajVec_Index_To_HitFit);

  // ana
  tree.SetBranchAddress("AnaVec_Matched_DetN", &pAnaVec_Matched_DetN);
  tree.SetBranchAddress("AnaVec_Matched_ResidU", &pAnaVec_Matched_ResdU);
  tree.SetBranchAddress("AnaVec_Matched_ResidV", &pAnaVec_Matched_ResdV);
}

size_t altel::TelEventTTreeReader::numEvents() const{
  return m_numEvents;
}

std::shared_ptr<altel::TelEvent> altel::TelEventTTreeReader::createTelEvent(size_t n){
  if(n>=m_numEvents){
    std::fprintf(stderr, "no more entries in ttree\n");
    return nullptr;
  }

  m_pTTree->GetEntry(n);
  // std::printf("GetEntry %d\n", n);
  // std::printf("rRunN %d, rEventN %d, rConfigN %d, rClock %d\n", rRunN, rEventN, rConfigN, rClock);
  // std::printf("measHit IDsize %d Usize %d Vsize %d\n", rHitMeasVec_DetN.size(), rHitMeasVec_U.size(), rHitMeasVec_V.size());
  // std::printf("fitHit IDsize %d Usize %d Vsize %d Xsize %d Ysize %d Zsize %d DXsize %d DYsize %d DZsize %d\n",
  //             rHitFitVec_DetN.size(), rHitFitVec_U.size(), rHitFitVec_V.size(),
  //             rHitFitVec_X.size(), rHitFitVec_Y.size(), rHitFitVec_Z.size(),
  //             rHitFitVec_DirX.size(), rHitFitVec_DirY.size(), rHitFitVec_DirZ.size());

  std::shared_ptr<altel::TelEvent> telEvent(new altel::TelEvent(rRunN, rEventN, rConfigN, rClock));


  auto it_rawMeasVec_DetN = rRawMeasVec_DetN.begin();
  auto it_rawMeasVec_U = rRawMeasVec_U.begin();
  auto it_rawMeasVec_V = rRawMeasVec_V.begin();
  auto it_rawMeasVec_Clk = rRawMeasVec_Clk.begin();
  auto it_rawMeasVec_DetN_end = rRawMeasVec_DetN.end();

  std::vector<altel::TelMeasRaw> measRaws;
  measRaws.reserve(rRawMeasVec_DetN.size());
  while(it_rawMeasVec_DetN !=it_rawMeasVec_DetN_end){
    measRaws.emplace_back(*it_rawMeasVec_U, *it_rawMeasVec_V, *it_rawMeasVec_DetN, *it_rawMeasVec_Clk);
    it_rawMeasVec_DetN++;
    it_rawMeasVec_U++;
    it_rawMeasVec_V++;
    it_rawMeasVec_Clk++;
  }

  // make index cluster
  auto it_numRawMeas_PerHitMeas = rHitMeasVec_NumRawMeas_PerHitMeas.begin();
  auto it_rawMeas_index = rHitMeasVec_Index_To_RawMeas.begin();

  auto it_numRawMeas_PerHitMeas_end = rHitMeasVec_NumRawMeas_PerHitMeas.end();
  auto it_rawMeas_index_end = rHitMeasVec_Index_To_RawMeas.end();

  std::vector<std::vector<altel::TelMeasRaw>> clusterCol;
  clusterCol.reserve(rHitMeasVec_DetN.size());
  while(it_numRawMeas_PerHitMeas != it_numRawMeas_PerHitMeas_end){
    std::vector<altel::TelMeasRaw> cluster;
    cluster.reserve(*it_numRawMeas_PerHitMeas);

    size_t numRawMeas_got = 0;
    while(it_rawMeas_index!=it_rawMeas_index_end && numRawMeas_got < *it_numRawMeas_PerHitMeas){
      int16_t rawMeasIndex = *it_rawMeas_index;
      if(rawMeasIndex!=int16_t(-1)){
        auto aMeasRaw = measRaws[rawMeasIndex];
        cluster.push_back(aMeasRaw);
        numRawMeas_got++;
      }
      it_rawMeas_index++;
    }
    clusterCol.push_back(std::move(cluster));
    it_numRawMeas_PerHitMeas++;
  }

  // std::cout<< "create cluster number "<<  clusterCol.size()<<std::endl;
  auto it_clusterCol = clusterCol.begin();

  size_t numMeasHit = rHitMeasVec_DetN.size();
  assert(rHitMeasVec_U.size() == numMeasHit && rHitMeasVec_V.size() == numMeasHit);
  auto it_measHitVec_detN = rHitMeasVec_DetN.begin();
  auto it_measHitVec_U = rHitMeasVec_U.begin();
  auto it_measHitVec_V = rHitMeasVec_V.begin();
  auto it_measHitVec_detN_end = rHitMeasVec_DetN.end();

  std::vector<std::shared_ptr<altel::TelMeasHit>> measHits;
  measHits.reserve(rHitMeasVec_DetN.size());
  while(it_measHitVec_detN !=it_measHitVec_detN_end){
    measHits.emplace_back(new altel::TelMeasHit(*it_measHitVec_detN,
                                                *it_measHitVec_U,
                                                *it_measHitVec_V,
                                                *it_clusterCol));
    // HitMeasVec_NumRawMeas_PerHitMeas;
    it_measHitVec_detN++;
    it_measHitVec_U++;
    it_measHitVec_V++;

    it_clusterCol++;
  }

  // std::cout<< "create measHits "<<  measHits.size()<<std::endl;

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

  // std::cout<<"loop fitHit"<<std::endl;
  while(it_fitHitVec_detN !=it_fitHitVec_detN_end){
    // std::cout<<"fitHitVec_detN "<< *it_fitHitVec_detN <<std::endl;
    int16_t originMeasHitIndex = *it_fitHitVec_OriginMeasHit_index;
    std::shared_ptr<altel::TelMeasHit> originMeasHit;
    // std::cout<< "originMeasHitIndex"<< originMeasHitIndex<<std::endl;
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


    int16_t matchedMeasHitIndex = *it_fitHitVec_MatchedMeasHit_index;
    std::shared_ptr<altel::TelMeasHit> matchedMeasHit;
    // std::cout<< "matchedMeasHitIndex"<< matchedMeasHitIndex<<std::endl;
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

  // std::cout<<"loop traj"<<std::endl;
  while(it_numFitHit_PerTraj != it_numFitHit_PerTraj_end){
    auto traj = std::make_shared<altel::TelTrajectory>();
    size_t numFitHits_got = 0;
    while(it_fitHit_index!=it_fitHit_index_end && numFitHits_got < *it_numFitHit_PerTraj){
      int16_t fitHitIndex = *it_fitHit_index;
      if(fitHitIndex!=int16_t(-1)){
        auto trajHit = trajHits[fitHitIndex];
        traj->trajHits().push_back(trajHit);
        numFitHits_got++;
      }
      it_fitHit_index++;
    }
    trajs.push_back(traj);
    it_numFitHit_PerTraj++;
  }

  telEvent->measHits() = std::move(measHits);
  telEvent->trajs() = std::move(trajs);

  return telEvent;
}
