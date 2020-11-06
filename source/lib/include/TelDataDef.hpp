#pragma once

#include <vector>
#include <memory>
#include <cstddef>

namespace TelActs{

  union TelMeasureRaw{ //64bits raw
    uint64_t mULL;
    int64_t  mLL;
    double   mD;

    int32_t  xyL[2];
    uint32_t xyUL[2];
    float    xyF[2];

    uint16_t xyztUS[4];
    int16_t  xyztS[4];

    char          vC[8];
    unsigned char vUC[8];
  };

  struct TelHitMeasure{
    uint32_t DN;    // detector id
    double PLs[2];  // local pos
    std::vector<TelMeasureRaw> Ms;// measures
  };

  struct TelHitFit{
    uint32_t DN;   // detector id
    double PGs[3]; // global pos
    double DGs[3]; // global dir

    double PLs[2]; // local pos
    double DLs[3]; // local dir
    std::shared_ptr<TelHitMeasure> HM; // measure hit
  };

  struct TelHit{
    uint32_t DN; // detecor id
    std::shared_ptr<TelHitFit> HF;    // fitted hit
    std::shared_ptr<TelHitMeasure> HM;// measure hit
  };

  struct TelTrajectory{
    uint32_t TN;                             // trajectories id
    std::vector<std::shared_ptr<TelHit>> Hs; // hit
  };

  struct TelEvent{
    uint32_t RN; // run id
    uint32_t EN; // event id
    std::vector<std::shared_ptr<TelTrajectory>> Ts; // trajectories
    std::vector<std::shared_ptr<TelHitMeasure>> HMs;// measure hit
  };
}


/*

{"E":
 {"EN":0,
  "RN":0,
  "DN":0,
  "Ms":[],

  "HMs":[
      {"HM":{"DN":0,"PLs":[],"Ms":[]}}
  ]

  "Ts":[
      {"T":
       {"TN":0,
        "Hs":[
            {"H":
             {"DN":0,
              "HF":{"DN":0, "MF": true, "PLs":[],"DLs":[],"PGs":[],"DGs":[]},
              "HM":{"DN":0, "PLs":[],"Ms":[]}
             }
            }
        ]
       }
  ]

 }
}

*/
