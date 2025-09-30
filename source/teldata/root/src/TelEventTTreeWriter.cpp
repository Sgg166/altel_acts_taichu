#include "TelEventTTreeWriter.hpp"
#include <TFile.h>
#include <TTree.h>

#include <iostream>

void altel::TelEventTTreeWriter::setTTree(TTree* pTTree){
  if(!pTTree){
    std::fprintf(stderr, "TTree is not yet set\n");
    throw;
  }
  m_pTTree = pTTree;
  TTree &tree = *m_pTTree;

  auto bRunN = tree.Branch("RunN", &rRunN);
  auto bEventN = tree.Branch("EveN", &rEventN);
  auto bConfigN = tree.Branch("DetN", &rConfigN);
  auto bClock = tree.Branch("ClkN", (ULong64_t*)&rClock);

  auto bNumTraj_PerEvent = tree.Branch("NumTrajs_PerEvent", &rNumTraj_PerEvent);
  auto bNumMeasHits_PerEvent = tree.Branch("NumMeasHits_PerEvent", &rNumMeasHits_PerEvent);

  auto bRawMeasVec_DetN = tree.Branch("MeasRawVec_DetN", &pRawMeasVec_DetN);
  auto bRawMeasVec_U = tree.Branch("MeasRawVec_U", &pRawMeasVec_U);
  auto bRawMeasVec_V = tree.Branch("MeasRawVec_V", &pRawMeasVec_V);
  auto bRawMeasVec_Clk = tree.Branch("MeasRawVec_Clk", &pRawMeasVec_Clk);

  auto bHitMeasVec_DetN = tree.Branch("MeasHitVec_DetN", &pHitMeasVec_DetN);
  auto bHitMeasVec_U = tree.Branch("MeasHitVec_U", &pHitMeasVec_U);
  auto bHitMeasVec_V = tree.Branch("MeasHitVec_V", &pHitMeasVec_V);
  auto bHitMeasVec_NumRawMeas_PerHitMeas =
    tree.Branch("MeasHitVec_NumMeasRaws_PerMeasHit", &pHitMeasVec_NumRawMeas_PerHitMeas);
  auto bHitMeasVec_Index_To_RawMeas =
    tree.Branch("MeasHitVec_MeasRaw_Index", &pHitMeasVec_Index_To_RawMeas);

  auto bHitFitVec_DetN = tree.Branch("TrajHitVec_DetN", &pHitFitVec_DetN);
  auto bHitFitVec_U = tree.Branch("TrajHitVec_U", &pHitFitVec_U);
  auto bHitFitVec_V = tree.Branch("TrajHitVec_V", &pHitFitVec_V);
  auto bHitFitVec_U_err = tree.Branch("TrajHitVec_U_err", &pHitFitVec_U_err);
  auto bHitFitVec_V_err = tree.Branch("TrajHitVec_V_err", &pHitFitVec_V_err);

  auto bHitFitVec_X = tree.Branch("TrajHitVec_X", &pHitFitVec_X);
  auto bHitFitVec_Y = tree.Branch("TrajHitVec_Y", &pHitFitVec_Y);
  auto bHitFitVec_Z = tree.Branch("TrajHitVec_Z", &pHitFitVec_Z);
  auto bHitFitVec_DirX = tree.Branch("TrajHitVec_DirX", &pHitFitVec_DirX);
  auto bHitFitVec_DirY = tree.Branch("TrajHitVec_DirY", &pHitFitVec_DirY);
  auto bHitFitVec_DirZ = tree.Branch("TrajHitVec_DirZ", &pHitFitVec_DirZ);
  auto bHitFitVec_Index_To_Origin_HitMeas =
    tree.Branch("TrajHitVec_OriginMeasHit_Index", &pHitFitVec_Index_To_Origin_HitMeas);
  auto bHitFitVec_Index_To_Matched_HitMeas =
    tree.Branch("TrajHitVec_MatchedMeasHit_Index", &pHitFitVec_Index_To_Matched_HitMeas);

  auto bTrajVec_NumHitFit_PerTraj = tree.Branch("TrajVec_NumTrajHits_PerTraj", &pTrajVec_NumHitFit_PerTraj);
  auto bTrajVec_NumHitMeas_Origin_PerTraj = tree.Branch("TrajVec_NumOriginMeasHits_PerTraj", &pTrajVec_NumHitMeas_Origin_PerTraj);
  auto bTrajVec_NumHitMeas_Matched_PerTraj = tree.Branch("TrajVec_NumMatchedMeasHits_PerTraj", &pTrajVec_NumHitMeas_Matched_PerTraj);

  auto bTrajVec_Index_To_HitFit =
    tree.Branch("TrajVec_TrajHit_Index", &pTrajVec_Index_To_HitFit);

  // ana
  auto bAnaVec_Matched_DetN = tree.Branch("AnaVec_Matched_DetN", &pAnaVec_Matched_DetN);
  auto bAnaVec_Matched_ResdU = tree.Branch("AnaVec_Matched_ResidU", &pAnaVec_Matched_ResdU);
  auto bAnaVec_Matched_ResdV = tree.Branch("AnaVec_Matched_ResidV", &pAnaVec_Matched_ResdV);
}


void altel::TelEventTTreeWriter::fillTelEvent(std::shared_ptr<altel::TelEvent> telEvent){

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
    
    rHitFitVec_U_err.clear();
    rHitFitVec_V_err.clear();
    
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
    std::map<altel::TelMeasRaw, int16_t> poolMapRawMeas;
    std::map<std::shared_ptr<altel::TelMeasHit>, int16_t> poolMapHitMeas;
    std::map<std::shared_ptr<altel::TelFitHit>, int16_t> poolMapHitFit;

    // for(auto &aRawMeas: telEvent->measRaws()){
    //   auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()) );
    //   if(inserted){
    //     // inserted and return true, when no exist
    //     rRawMeasVec_U.push_back(aRawMeas.u());
    //     rRawMeasVec_V.push_back(aRawMeas.v());
    //     rRawMeasVec_DetN.push_back(aRawMeas.detN());
    //     rRawMeasVec_Clk.push_back(aRawMeas.clkN());
    //   }
    // }

    for(auto &aHitMeas: telEvent->measHits()){
      auto [it, inserted] = poolMapHitMeas.emplace( aHitMeas, int16_t(rHitMeasVec_DetN.size()) );
      if(inserted){
        // inserted and return true, when no exist
        rHitMeasVec_DetN.push_back(aHitMeas->detN());
        rHitMeasVec_U.push_back(aHitMeas->u());
        rHitMeasVec_V.push_back(aHitMeas->v());

        rHitMeasVec_Index_To_RawMeas.push_back(-1);//: -1 {M-N} -1 {M-N}
        for(auto &aRawMeas: aHitMeas->measRaws()){
          auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()));
          if(inserted){
            rRawMeasVec_U.push_back(aRawMeas.u());
            rRawMeasVec_V.push_back(aRawMeas.v());
            rRawMeasVec_DetN.push_back(aRawMeas.detN());
            rRawMeasVec_Clk.push_back(aRawMeas.clkN());
          }
          rHitMeasVec_Index_To_RawMeas.push_back(it->second);
        }
        rHitMeasVec_NumRawMeas_PerHitMeas.push_back(int16_t(aHitMeas->measRaws().size()));
      }
    }

    rNumTraj_PerEvent = 0;
    for(auto &aTraj: telEvent->trajs()){
      size_t fittedHitNum = aTraj->numOriginMeasHit();
      // if(fittedHitNum<5){
      //   continue;
      // }

      rNumTraj_PerEvent ++;
      rTrajVec_Index_To_HitFit.push_back(-1);//traj: -1 {M-N} -1 {M-N}

      int16_t aNumHitFit_PerTraj = 0;
      int16_t aNumHitMeas_Origin_PerTraj = 0;
      int16_t aNumHitMeas_Matched_PerTraj = 0;

      for(auto &aHit: aTraj->trajHits()){
        aNumHitFit_PerTraj++;
        auto aHitFit = aHit->fitHit();
        auto [it, inserted] = poolMapHitFit.emplace( aHitFit, int16_t(rHitFitVec_DetN.size()) );
        if(!inserted){
          //TODO: handle this case for ambiguility case
          std::fprintf(stderr, "a shared hitfit by multi-traj, WRONG\n");
          throw;
        }

        rTrajVec_Index_To_HitFit.push_back(it->second);
        rHitFitVec_TrajN.push_back(int16_t(rTrajVec_NumHitFit_PerTraj.size()));
        rHitFitVec_DetN.push_back(aHitFit->detN());
        rHitFitVec_U.push_back(aHitFit->u());
        rHitFitVec_V.push_back(aHitFit->v());

        rHitFitVec_U_err.push_back(aHitFit->u_err());
        rHitFitVec_V_err.push_back(aHitFit->v_err());
        
        rHitFitVec_X.push_back(aHitFit->x());
        rHitFitVec_Y.push_back(aHitFit->y());
        rHitFitVec_Z.push_back(aHitFit->z());
        rHitFitVec_DirX.push_back(aHitFit->dx());
        rHitFitVec_DirY.push_back(aHitFit->dy());
        rHitFitVec_DirZ.push_back(aHitFit->dz());

        int16_t indexHitMeasOri = -1;
        auto aHitMeasOri = aHitFit->OM;
        if(aHitMeasOri){
          aNumHitMeas_Origin_PerTraj ++;
          auto [it, inserted] = poolMapHitMeas.emplace( aHitMeasOri, int16_t(rHitMeasVec_DetN.size()));
          if(inserted){
            rHitMeasVec_DetN.push_back(aHitMeasOri->detN());
            rHitMeasVec_U.push_back(aHitMeasOri->u());
            rHitMeasVec_V.push_back(aHitMeasOri->v());

            rHitMeasVec_Index_To_RawMeas.push_back(-1);//: -1 {M-N} -1 {M-N}
            for(auto &aRawMeas: aHitMeasOri->measRaws()){
              auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()));
              if(inserted){
                rRawMeasVec_U.push_back(aRawMeas.u());
                rRawMeasVec_V.push_back(aRawMeas.v());
                rRawMeasVec_DetN.push_back(aRawMeas.detN());
                rRawMeasVec_Clk.push_back(aRawMeas.clkN());
              }
              rHitMeasVec_Index_To_RawMeas.push_back(it->second);
            }
            rHitMeasVec_NumRawMeas_PerHitMeas.push_back(int16_t(aHitMeasOri->measRaws().size()));
          }
          indexHitMeasOri = it->second;
        }
        rHitFitVec_Index_To_Origin_HitMeas.push_back(indexHitMeasOri);

        int16_t  indexHitMeasMatched = -1;
        auto aHitMeasMatched = aHit->MM;
        if(aHitMeasMatched){
          aNumHitMeas_Matched_PerTraj++;
          auto [it, inserted] = poolMapHitMeas.emplace( aHitMeasMatched, int16_t(rHitMeasVec_DetN.size()));
          if(inserted){
            rHitMeasVec_DetN.push_back(aHitMeasMatched->detN());
            rHitMeasVec_U.push_back(aHitMeasMatched->u());
            rHitMeasVec_V.push_back(aHitMeasMatched->v());

            rHitMeasVec_Index_To_RawMeas.push_back(-1);//: -1 {M-N} -1 {M-N}
            for(auto &aRawMeas: aHitMeasMatched->measRaws()){
              auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()));
              if(inserted){
                rRawMeasVec_U.push_back(aRawMeas.u());
                rRawMeasVec_V.push_back(aRawMeas.v());
                rRawMeasVec_DetN.push_back(aRawMeas.detN());
                rRawMeasVec_Clk.push_back(aRawMeas.clkN());
              }
              rHitMeasVec_Index_To_RawMeas.push_back(it->second);
            }
            rHitMeasVec_NumRawMeas_PerHitMeas.push_back(int16_t(aHitMeasMatched->measRaws().size()));
            //

          }
          indexHitMeasMatched = it->second;
        }
        rHitFitVec_Index_To_Matched_HitMeas.push_back(indexHitMeasMatched);

        //ana matched
        if(aHitMeasMatched){
          rAnaVec_Matched_DetN.push_back(aHitMeasMatched->detN());
          rAnaVec_Matched_ResdU.push_back(aHitMeasMatched->u() - aHitFit->u());
          rAnaVec_Matched_ResdV.push_back(aHitMeasMatched->v() - aHitFit->v());
        }
      }
      rTrajVec_NumHitFit_PerTraj.push_back(aNumHitFit_PerTraj);
      rTrajVec_NumHitMeas_Origin_PerTraj.push_back(aNumHitMeas_Origin_PerTraj);
      rTrajVec_NumHitMeas_Matched_PerTraj.push_back(aNumHitMeas_Matched_PerTraj);
    }

    rRunN = telEvent->runN();
    rEventN = telEvent->eveN();
    rConfigN = telEvent->detN();
    rClock = telEvent->clkN();
    rNumMeasHits_PerEvent = poolMapHitMeas.size();

    m_pTTree->Fill();
}
