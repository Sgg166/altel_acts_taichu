#include "TelEventTTreeWriter.hpp"
#include <TFile.h>
#include <TTree.h>

#include <iostream>

void TelActs::TelEventTTreeWriter::createBranch(){
  if(!pTree){
    std::fprintf(stderr, "TTree is not yet set\n");
    throw;
  }
  TTree &tree = *pTree;

  auto bRunN = tree.Branch("RunN", &rRunN);
  auto bEventN = tree.Branch("EventN", &rEventN);
  auto bConfigN = tree.Branch("ConfigN", &rConfigN);
  auto bNumTraj_PerEvent = tree.Branch("NumTraj_PerEvent", &rNumTraj_PerEvent);

  auto bRawMeasVec_DetN = tree.Branch("RawMeasVec_DetN", &pRawMeasVec_DetN);
  auto bRawMeasVec_U = tree.Branch("RawMeasVec_U", &pRawMeasVec_U);
  auto bRawMeasVec_V = tree.Branch("RawMeasVec_V", &pRawMeasVec_V);
  auto bRawMeasVec_Clk = tree.Branch("RawMeasVec_Clk", &pRawMeasVec_Clk);

  auto bHitMeasVec_DetN = tree.Branch("HitMeasVec_DetN", &pHitMeasVec_DetN);
  auto bHitMeasVec_U = tree.Branch("HitMeasVec_U", &pHitMeasVec_U);
  auto bHitMeasVec_V = tree.Branch("HitMeasVec_V", &pHitMeasVec_V);
  auto bHitMeasVec_NumRawMeas_PerHitMeas =
    tree.Branch("HitMeasVec_NumRawMeas_PerHitMeas", &pHitMeasVec_NumRawMeas_PerHitMeas);
  auto bHitMeasVec_Index_To_RawMeas =
    tree.Branch("HitMeasVec_Index_To_RawMeas", &pHitMeasVec_Index_To_RawMeas);

  auto bHitFitVec_DetN = tree.Branch("HitFitVec_DetN", &pHitFitVec_DetN);
  auto bHitFitVec_U = tree.Branch("HitFitVec_U", &pHitFitVec_U);
  auto bHitFitVec_V = tree.Branch("HitFitVec_V", &pHitFitVec_V);
  auto bHitFitVec_X = tree.Branch("HitFitVec_X", &pHitFitVec_X);
  auto bHitFitVec_Y = tree.Branch("HitFitVec_Y", &pHitFitVec_Y);
  auto bHitFitVec_Z = tree.Branch("HitFitVec_Z", &pHitFitVec_Z);
  auto bHitFitVec_DirX = tree.Branch("HitFitVec_DirX", &pHitFitVec_DirX);
  auto bHitFitVec_DirY = tree.Branch("HitFitVec_DirY", &pHitFitVec_DirY);
  auto bHitFitVec_DirZ = tree.Branch("HitFitVec_DirZ", &pHitFitVec_DirZ);
  auto bHitFitVec_Index_To_Origin_HitMeas =
    tree.Branch("HitFitVec_Index_To_Origin_HitMeas", &pHitFitVec_Index_To_Origin_HitMeas);
  auto bHitFitVec_Index_To_Matched_HitMeas =
    tree.Branch("HitFitVec_Index_To_Matched_HitMeas", &pHitFitVec_Index_To_Matched_HitMeas);

  auto bTrajVec_NumHitFit_PerTraj = tree.Branch("TrajVec_NumHitFit_PerTraj", &pTrajVec_NumHitFit_PerTraj);
  auto bTrajVec_NumHitMeas_Origin_PerTraj = tree.Branch("TrajVec_NumHitMeas_Origin_PerTraj", &pTrajVec_NumHitMeas_Origin_PerTraj);
  auto bTrajVec_NumHitMeas_Matched_PerTraj = tree.Branch("TrajVec_NumHitMeas_Matched_PerTraj", &pTrajVec_NumHitMeas_Matched_PerTraj);

  auto bTrajVec_Index_To_HitFit =
    tree.Branch("TrajVec_Index_To_HitFit", &pTrajVec_Index_To_HitFit);

  // ana
  auto bAnaVec_Matched_DetN = tree.Branch("AnaVec_Matched_DetN", &pAnaVec_Matched_DetN);
  auto bAnaVec_Matched_ResdU = tree.Branch("AnaVec_Matched_ResdU", &pAnaVec_Matched_ResdU);
  auto bAnaVec_Matched_ResdV = tree.Branch("AnaVec_Matched_ResdV", &pAnaVec_Matched_ResdV);
}


void TelActs::TelEventTTreeWriter::fill(std::shared_ptr<TelActs::TelEvent> telEvent){

///////////////branch vector clear////////////////////////
    ///rawMeas
    rRawMeasVec_DetN.clear();
    rRawMeasVec_U.clear();
    rRawMeasVec_V.clear();
    rRawMeasVec_Clk.clear();
    //hitMeas
    rHitMeasVec_DetN.clear();
    rHitMeasVec_U.clear();
    rHitMeasVec_V.clear();
    rHitMeasVec_NumRawMeas_PerHitMeas.clear();
    rHitMeasVec_Index_To_RawMeas.clear();

    //hitFit
    rHitFitVec_TrajN.clear();
    rHitFitVec_DetN.clear();
    rHitFitVec_U.clear();
    rHitFitVec_V.clear();
    rHitFitVec_X.clear();
    rHitFitVec_Y.clear();
    rHitFitVec_Z.clear();
    rHitFitVec_DirX.clear();
    rHitFitVec_DirY.clear();
    rHitFitVec_DirZ.clear();
    rHitFitVec_Index_To_Origin_HitMeas.clear();
    rHitFitVec_Index_To_Matched_HitMeas.clear();

    //traj
    rTrajVec_NumHitFit_PerTraj.clear();
    rTrajVec_NumHitMeas_Origin_PerTraj.clear();
    rTrajVec_NumHitMeas_Matched_PerTraj.clear();
    rTrajVec_Index_To_HitFit.clear();


    //matched analysis
    rAnaVec_Matched_DetN.clear();
    rAnaVec_Matched_ResdU.clear();
    rAnaVec_Matched_ResdV.clear();
///////////////////////////////////////////////


////////////////////////////////////////fill tree////////////////////
    std::map<std::shared_ptr<TelActs::TelRawMeasure>, int16_t> poolMapRawMeas;
    std::map<std::shared_ptr<TelActs::TelHitMeasure>, int16_t> poolMapHitMeas;
    std::map<std::shared_ptr<TelActs::TelHitFit>, int16_t> poolMapHitFit;

    for(auto &aRawMeas: telEvent->Ms){
      auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()) );
      if(inserted){
        // inserted and return true, when no exist
        rRawMeasVec_U.push_back(aRawMeas->uvdcS[0]);
        rRawMeasVec_V.push_back(aRawMeas->uvdcS[1]);
        rRawMeasVec_DetN.push_back(aRawMeas->uvdcS[2]);
        rRawMeasVec_Clk.push_back(aRawMeas->uvdcS[3]);
      }
    }

    for(auto &aHitMeas: telEvent->HMs){
      auto [it, inserted] = poolMapHitMeas.emplace( aHitMeas, int16_t(rHitMeasVec_DetN.size()) );
      if(inserted){
        // inserted and return true, when no exist
        rHitMeasVec_DetN.push_back(aHitMeas->DN);
        rHitMeasVec_U.push_back(aHitMeas->PLs[0]);
        rHitMeasVec_V.push_back(aHitMeas->PLs[1]);

        rHitMeasVec_Index_To_RawMeas.push_back(-1);//: -1 {M-N} -1 {M-N}
        for(auto &aRawMeas: aHitMeas->Ms){
          auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()));
          if(inserted){
            rRawMeasVec_U.push_back(aRawMeas->uvdcS[0]);
            rRawMeasVec_V.push_back(aRawMeas->uvdcS[1]);
            rRawMeasVec_DetN.push_back(aRawMeas->uvdcS[2]);
            rRawMeasVec_Clk.push_back(aRawMeas->uvdcS[3]);
          }
          rHitMeasVec_Index_To_RawMeas.push_back(it->second);
        }
        rHitMeasVec_NumRawMeas_PerHitMeas.push_back(int16_t(aHitMeas->Ms.size()));
      }
    }

    rNumTraj_PerEvent = 0;
    for(auto &aTraj: telEvent->Ts){
      size_t fittedHitNum = aTraj->numberHitFitByMeas();
      if(fittedHitNum<5){
        continue;
      }

      rNumTraj_PerEvent ++;
      rTrajVec_Index_To_HitFit.push_back(-1);//traj: -1 {M-N} -1 {M-N}

      int16_t aNumHitFit_PerTraj = 0;
      int16_t aNumHitMeas_Origin_PerTraj = 0;
      int16_t aNumHitMeas_Matched_PerTraj = 0;

      for(auto &aHit: aTraj->Hs){
        aNumHitFit_PerTraj++;
        auto aHitFit = aHit->HF;
        auto [it, inserted] = poolMapHitFit.emplace( aHitFit, int16_t(rHitFitVec_DetN.size()) );
        if(!inserted){
          //TODO: handle this case for ambiguility case
          std::fprintf(stderr, "a shared hitfit by multi-traj, WRONG\n");
          throw;
        }

        rTrajVec_Index_To_HitFit.push_back(it->second);
        rHitFitVec_TrajN.push_back(int16_t(rTrajVec_NumHitFit_PerTraj.size()));
        rHitFitVec_DetN.push_back(aHitFit->DN);
        rHitFitVec_U.push_back(aHitFit->PLs[0]);
        rHitFitVec_V.push_back(aHitFit->PLs[1]);
        rHitFitVec_X.push_back(aHitFit->PGs[0]);
        rHitFitVec_Y.push_back(aHitFit->PGs[1]);
        rHitFitVec_Z.push_back(aHitFit->PGs[2]);
        rHitFitVec_DirX.push_back(aHitFit->DGs[0]);
        rHitFitVec_DirY.push_back(aHitFit->DGs[1]);
        rHitFitVec_DirZ.push_back(aHitFit->DGs[2]);

        int16_t indexHitMeasOri = -1;
        auto aHitMeasOri = aHitFit->HM;
        if(aHitMeasOri){
          aNumHitMeas_Origin_PerTraj ++;
          auto [it, inserted] = poolMapHitMeas.emplace( aHitMeasOri, int16_t(rHitMeasVec_DetN.size()));
          if(inserted){
            rHitMeasVec_DetN.push_back(aHitMeasOri->DN);
            rHitMeasVec_U.push_back(aHitMeasOri->PLs[0]);
            rHitMeasVec_V.push_back(aHitMeasOri->PLs[1]);

            rHitMeasVec_Index_To_RawMeas.push_back(-1);//: -1 {M-N} -1 {M-N}
            for(auto &aRawMeas: aHitMeasOri->Ms){
              auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()));
              if(inserted){
                rRawMeasVec_U.push_back(aRawMeas->uvdcS[0]);
                rRawMeasVec_V.push_back(aRawMeas->uvdcS[1]);
                rRawMeasVec_DetN.push_back(aRawMeas->uvdcS[2]);
                rRawMeasVec_Clk.push_back(aRawMeas->uvdcS[3]);
              }
              rHitMeasVec_Index_To_RawMeas.push_back(it->second);
            }
            rHitMeasVec_NumRawMeas_PerHitMeas.push_back(int16_t(aHitMeasOri->Ms.size()));
          }
          indexHitMeasOri = it->second;
        }
        rHitFitVec_Index_To_Origin_HitMeas.push_back(indexHitMeasOri);

        int16_t  indexHitMeasMatched = -1;
        auto aHitMeasMatched = aHit->HM;
        if(aHitMeasMatched){
          aNumHitMeas_Matched_PerTraj++;
          auto [it, inserted] = poolMapHitMeas.emplace( aHitMeasMatched, int16_t(rHitMeasVec_DetN.size()));
          if(inserted){
            rHitMeasVec_DetN.push_back(aHitMeasMatched->DN);
            rHitMeasVec_U.push_back(aHitMeasMatched->PLs[0]);
            rHitMeasVec_V.push_back(aHitMeasMatched->PLs[1]);

            rHitMeasVec_Index_To_RawMeas.push_back(-1);//: -1 {M-N} -1 {M-N}
            for(auto &aRawMeas: aHitMeasMatched->Ms){
              auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()));
              if(inserted){
                rRawMeasVec_U.push_back(aRawMeas->uvdcS[0]);
                rRawMeasVec_V.push_back(aRawMeas->uvdcS[1]);
                rRawMeasVec_DetN.push_back(aRawMeas->uvdcS[2]);
                rRawMeasVec_Clk.push_back(aRawMeas->uvdcS[3]);
              }
              rHitMeasVec_Index_To_RawMeas.push_back(it->second);
            }
            rHitMeasVec_NumRawMeas_PerHitMeas.push_back(int16_t(aHitMeasMatched->Ms.size()));
            //

          }
          indexHitMeasMatched = it->second;
        }
        rHitFitVec_Index_To_Matched_HitMeas.push_back(indexHitMeasMatched);

        //ana matched
        if(aHitMeasMatched){
          rAnaVec_Matched_DetN.push_back(aHitMeasMatched->DN);
          rAnaVec_Matched_ResdU.push_back(aHitMeasMatched->PLs[0] - aHitFit->PLs[0]);
          rAnaVec_Matched_ResdV.push_back(aHitMeasMatched->PLs[1] - aHitFit->PLs[1]);
        }
      }
      rTrajVec_NumHitFit_PerTraj.push_back(aNumHitFit_PerTraj);
      rTrajVec_NumHitMeas_Origin_PerTraj.push_back(aNumHitMeas_Origin_PerTraj);
      rTrajVec_NumHitMeas_Matched_PerTraj.push_back(aNumHitMeas_Matched_PerTraj);
    }

    rRunN = telEvent->RN;
    rEventN = telEvent->EN;
    rConfigN = telEvent->DN;

    // if(!pTree){
    //   std::fprintf(stderr, "TTree is not yet set\n");
    //   throw;
    // }
    // TTree &tree = *pTree;
    pTree->Fill();
    // std::cout<<"reach ptree->fill"<<std::endl;
}
