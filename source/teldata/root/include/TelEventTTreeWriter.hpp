#pragma once

#include "TelEvent.hpp"


//TODO TTree::MakeClass  https://root.cern.ch/doc/master/classTTree.html#ac4ceaf4ae0b87412acf94093043cc2de

class TFile;
class TTree;

namespace altel{
  struct TelEventTTreeWriter{
    uint32_t rRunN;
    uint32_t rEventN;
    int16_t rConfigN;
    uint64_t rClock; // timestamp / trigger id;

    int16_t rNumTraj_PerEvent;

    ///rawMeas
    std::vector<int16_t> rRawMeasVec_DetN, *pRawMeasVec_DetN = &rRawMeasVec_DetN;
    std::vector<int16_t> rRawMeasVec_U, *pRawMeasVec_U = &rRawMeasVec_U;
    std::vector<int16_t> rRawMeasVec_V, *pRawMeasVec_V = &rRawMeasVec_V;
    std::vector<int16_t> rRawMeasVec_Clk, *pRawMeasVec_Clk = &rRawMeasVec_Clk;

    //hitMeas
    std::vector<int16_t> rHitMeasVec_DetN, *pHitMeasVec_DetN = &rHitMeasVec_DetN;
    std::vector<double>  rHitMeasVec_U, *pHitMeasVec_U = &rHitMeasVec_U;
    std::vector<double>  rHitMeasVec_V, *pHitMeasVec_V = &rHitMeasVec_V;
    std::vector<int16_t> rHitMeasVec_NumRawMeas_PerHitMeas,
      *pHitMeasVec_NumRawMeas_PerHitMeas = &rHitMeasVec_NumRawMeas_PerHitMeas;
    std::vector<int16_t> rHitMeasVec_Index_To_RawMeas,
      *pHitMeasVec_Index_To_RawMeas = &rHitMeasVec_Index_To_RawMeas;

    //hitFit
    std::vector<int16_t> rHitFitVec_TrajN, *pHitFitVec_TrajN = &rHitFitVec_TrajN;
    std::vector<int16_t> rHitFitVec_DetN, *pHitFitVec_DetN = &rHitFitVec_DetN;
    std::vector<double> rHitFitVec_U, *pHitFitVec_U = &rHitFitVec_U;
    std::vector<double> rHitFitVec_V, *pHitFitVec_V = &rHitFitVec_V;
    std::vector<double> rHitFitVec_X, *pHitFitVec_X = &rHitFitVec_X;
    std::vector<double> rHitFitVec_Y, *pHitFitVec_Y = &rHitFitVec_Y;
    std::vector<double> rHitFitVec_Z, *pHitFitVec_Z = &rHitFitVec_Z;
    std::vector<double> rHitFitVec_DirX, *pHitFitVec_DirX = &rHitFitVec_DirX;
    std::vector<double> rHitFitVec_DirY, *pHitFitVec_DirY = &rHitFitVec_DirY;
    std::vector<double> rHitFitVec_DirZ, *pHitFitVec_DirZ = &rHitFitVec_DirZ;

    std::vector<int16_t> rHitFitVec_Index_To_Origin_HitMeas,
      *pHitFitVec_Index_To_Origin_HitMeas = &rHitFitVec_Index_To_Origin_HitMeas;

    std::vector<int16_t> rHitFitVec_Index_To_Matched_HitMeas,
      *pHitFitVec_Index_To_Matched_HitMeas = &rHitFitVec_Index_To_Matched_HitMeas;

    //traj
    std::vector<int16_t> rTrajVec_NumHitFit_PerTraj,
      *pTrajVec_NumHitFit_PerTraj = &rTrajVec_NumHitFit_PerTraj;
    std::vector<int16_t> rTrajVec_NumHitMeas_Origin_PerTraj,
      *pTrajVec_NumHitMeas_Origin_PerTraj = &rTrajVec_NumHitMeas_Origin_PerTraj;
    std::vector<int16_t> rTrajVec_NumHitMeas_Matched_PerTraj,
      *pTrajVec_NumHitMeas_Matched_PerTraj = &rTrajVec_NumHitMeas_Matched_PerTraj;
    std::vector<int16_t> rTrajVec_Index_To_HitFit,
      *pTrajVec_Index_To_HitFit = &rTrajVec_Index_To_HitFit;

    //ana
    std::vector<int16_t> rAnaVec_Matched_DetN, *pAnaVec_Matched_DetN = &rAnaVec_Matched_DetN;
    std::vector<double> rAnaVec_Matched_ResdU, *pAnaVec_Matched_ResdU = &rAnaVec_Matched_ResdU;
    std::vector<double> rAnaVec_Matched_ResdV, *pAnaVec_Matched_ResdV = &rAnaVec_Matched_ResdV;

    TTree* pTree;
    void createBranch();
    void fill(std::shared_ptr<altel::TelEvent> ev);
  };
}
