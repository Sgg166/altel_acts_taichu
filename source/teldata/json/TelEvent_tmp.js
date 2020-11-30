
{"EV":
 {"RN":runN64,
  "EN":eventN64,
  "DN":detSetupN16,
  "CK":clk64,
  "MRs":[[u16, v16, detN16, clk16]],
  "MHs":[
      {"MH":{"DN":detN16,"PLs":[u, v],"MRs":[[u16, v16, detN16, clkN16]]}}
  ],

  "TJs":[
      {"TJ":
       {"TN":trajN64,
        "THs":[
            {"TH":
             {"DN":detN16,
              "FH":{"DN":detN16,
                    "PLs":[u, v],"DLs":[du, dv, dw],
                    "PGs":[x, y, z],"DGs":[dx, dy, dz]
                   },
              "OM":{"MHi":MHs_index},
              "MM":{"MHi":MHs_index}
             }
            }
        ]
       }
      }
  ]
 }
}

