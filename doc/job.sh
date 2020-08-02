# 1) Align to get alignment geometry
./AlignTelescopeGeometry -f data/data_2GeV.json -e 2 -out geometry/geometry_2GeV_1.json
./AlignTelescopeGeometry -f data/data_3GeV.json -e 3 -out geometry/geometry_3GeV_1.json
./AlignTelescopeGeometry -f data/data_4GeV.json -e 4 -out geometry/geometry_4GeV_1.json
./AlignTelescopeGeometry -f data/data_5GeV.json -e 5 -out geometry/geometry_5GeV_1.json
./AlignTelescopeGeometry -f data/data_6GeV.json -e 6 -out geometry/geometry_6GeV_1.json

# (Optional) Perform Kalman fit with aligned geometry (only good tracks)
./FitTelescopeTracks -f data/data_2GeV.json  -e 2 -g geometry/geometry_2GeV_1.json  -o performance_fit_2GeV
./FitTelescopeTracks -f data/data_3GeV.json  -e 3 -g geometry/geometry_3GeV_1.json  -o performance_fit_3GeV
./FitTelescopeTracks -f data/data_4GeV.json  -e 4 -g geometry/geometry_4GeV_1.json  -o performance_fit_4GeV
./FitTelescopeTracks -f data/data_5GeV.json  -e 5 -g geometry/geometry_5GeV_1.json  -o performance_fit_5GeV
./FitTelescopeTracks -f data/data_6GeV.json  -e 6 -g geometry/geometry_6GeV_1.json  -o performance_fit_6GeV

# 2) Perform CKF track fitting/finding with aligned geometry (for all tracks)
./RecTelescopeCKFTracks -f data/data_2GeV.json  -e 2 -g geometry/geometry_2GeV_1.json -o performance_ckf_2GeV

# 3) Visualization of one event (measurements on found track(s) and tracking geometry)
meshlab Tracker_sensitives_l1.obj  Tracker_sensitives_l3.obj Tracker_sensitives_l5.obj Tracker_sensitives_l7.obj Tracker_sensitives_l9.obj Tracker_sensitives_l11.obj performance_ckf_2GeV/event000002103-TelescopeTrack.obj

# Print event-level raw hits
./PrintEventJson -f data/data_2GeV.json -e 200

