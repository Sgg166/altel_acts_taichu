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


  using TelHitMeas = TelHitMeasure;
  using TelTraj = TelTrajectory;

  using TelRawMeas = TelRawMeasure;
  using TelMeasRaw = TelRawMeasure;
  using TelMeasHit = TelHitMeasure;
  using TelFitHit = TelHitFit;
  using TelTrajHit = TelHit;

  struct TelEvent{
    uint64_t RN{0}; // run id
    uint64_t EN{0}; // event id
    uint64_t DN{0}; // detector/setup id
    std::vector<std::shared_ptr<TelRawMeasure>> Ms;// measures raw
    std::vector<std::shared_ptr<TelHitMeasure>> HMs;// measure hit
    std::vector<std::shared_ptr<TelTrajectory>> Ts; // trajectories

  };

}
