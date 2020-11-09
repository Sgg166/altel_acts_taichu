#pragma once

#include <vector>
#include <memory>
#include <cstddef>

namespace TelActs{

  union TelRawMeasure{ //64bits raw
    uint64_t mULL;
    int64_t  mLL;
    double   mD;

    int32_t  uvL[2];
    uint32_t uvUL[2];
    float    uvF[2];

    uint16_t uvdcUS[4];
    int16_t  uvdcS[4];

    char          vC[8];
    unsigned char vUC[8];
  };

  struct TelHitMeasure{
    uint64_t DN{0};    // detector id
    double PLs[2]{0.0, 0.0};  // local pos
    std::vector<std::shared_ptr<TelRawMeasure>> Ms;// measures
  };

  struct TelHitFit{
    uint64_t DN{0};   // detector id
    double PLs[2]{0.0, 0.0}; // local pos
    double DLs[3]{0.0, 0.0, 0.0}; // local dir

    double PGs[3]{0.0, 0.0, 0.0}; // global pos
    double DGs[3]{0.0, 0.0, 0.0}; // global dir

    std::shared_ptr<TelHitMeasure> HM; // measure hit

    bool isFittedFromMeasure() const{
      if(HM){
        return true;
      }
      else{
        return false;
      }
    }
  };

  struct TelHit{
    uint64_t DN{0}; // detector id
    std::shared_ptr<TelHitFit> HF;    // fitted hit
    std::shared_ptr<TelHitMeasure> HM;// measure hit

    bool hasHitFit() const{
      if(HF){
        return true;
      }
      else{
        return false;
      }
    }

    bool isFittedFromMeasure() const{
      if( HF && HF->isFittedFromMeasure()){
        return true;
      }
      else{
        return false;
      }
    }
  };

  struct TelTrajectory{
    uint64_t TN{0};                          // trajectories id
    std::vector<std::shared_ptr<TelHit>> Hs; // hit

    std::shared_ptr<TelHit> hit(size_t id){
      for(const auto& h: Hs){
        if (h && h->DN == id){
          return h;
        }
      }
      return nullptr;
    }

    size_t numberHitFitByMeas()const{
      size_t n = 0;
      for(const auto& h: Hs){
        if ( h && h->HF && h->HF->HM){
          n++;
        }
      }
      return n;
    }

  };

  struct TelEvent{
    uint64_t RN{0}; // run id
    uint64_t EN{0}; // event id
    uint64_t DN{0}; // detector/setup id
    std::vector<std::shared_ptr<TelRawMeasure>> Ms;// measures
    std::vector<std::shared_ptr<TelHitMeasure>> HMs;// measure hit
    std::vector<std::shared_ptr<TelTrajectory>> Ts; // trajectories

  };
}


/*

{"E":
 {"RN":0,
  "EN":0,
  "DN":0,

  "Ms":[],

  "HMs":[
      {"HM":{"DN":0,"PLs":[<double>, <double>],"Ms":[<int64>]}}
  ],

  "Ts":[
      {"T":
       {"TN":0,
        "Hs":[
            {"H":
             {"DN":0,
              "HF":{"DN":0, "MF": boolean,
                    "PLs":[<double>, <double>],"DLs":[<double>, <double>, <double>],
                    "PGs":[<double>, <double>, <double>],"DGs":[<double>, <double>, <double>]},
              "HM":{"DN":0, "PLs":[<double>, <double>],"Ms":[<int64>]}
             }
            }
        ]
       }
  ]

 }
}

*/


/* root branch

RunN
EventN
ConfigN

RawMeasVec_U
RawMeasVec_V
RawMeasVec_DetN

HitMeasVec_U
HitMeasVec_V
HitMeasVec_DetN
HitMeasVec_Index_to_RawMeas

HitFitVec_U
HitFitVec_V
HitFitVec_X
HitFitVec_Y
HitFitVec_Z
HitFitVec_DirX
HitFitVec_DirY
HitFitVec_DirZ
HitFitVec_DetN
HitFitVec_Bool_Has_Origin_HitMeas
HitFitVec_Bool_Has_Matched_HitMeas
HitFitVec_Index_to_Origin_HitMeas
HitFitVec_Index_to_Matched_HitMeas


TrajVec_IndexVec_to_HitFit

*/


