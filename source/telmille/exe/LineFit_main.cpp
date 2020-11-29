#include "getopt.h"
#include "mysystem.hh"
#include "myrapidjson.h"

#include <iostream>
#include <algorithm>
#include <regex>

#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH2.h>
#include <TImage.h>
#include <TCanvas.h>

template<typename Vector3>
std::pair<Vector3, Vector3> best_plane_from_points(const std::vector<Vector3> & c)
{
  // copy coordinates to  matrix in Eigen format
  size_t num_atoms = c.size();
  Eigen::Matrix< typename Vector3::Scalar, Eigen::Dynamic, Eigen::Dynamic > coord(3, num_atoms);
  for (size_t i = 0; i < num_atoms; ++i){
    coord.col(i) = c[i];
  }

  // calculate centroid
  Vector3 centroid(coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());

  // subtract centroid
  coord.row(0).array() -= centroid(0);
  coord.row(1).array() -= centroid(1);
  coord.row(2).array() -= centroid(2);

  // we only need the left-singular matrix here
  //  http://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
  auto svd = coord.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  // Vector3 plane_normal = svd.matrixU().rightCols<1>();
  Vector3 plane_normal = svd.matrixU().rightCols(1);
  return std::make_pair(centroid, plane_normal);
}


template<typename Vector3>
std::pair < Vector3, Vector3 > best_line_from_points(const std::vector<Vector3> & c)
{
  // copy coordinates to  matrix in Eigen format
  size_t num_atoms = c.size();
  Eigen::Matrix< typename Vector3::Scalar, Eigen::Dynamic, Eigen::Dynamic > centers(num_atoms, 3);
  for (size_t i = 0; i < num_atoms; ++i){
    centers.row(i) = c[i];
  }

  Vector3 origin = centers.colwise().mean();
  Eigen::MatrixXd centered = centers.rowwise() - origin.transpose();
  Eigen::MatrixXd cov = centered.adjoint() * centered;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
  Vector3 axis = eig.eigenvectors().col(2).normalized();

  return std::make_pair(origin, axis);
}


void FitTrack(unsigned int nMeasures,
                          const std::vector<double>& xPosMeasure,
                          const std::vector<double>& yPosMeasure,
                          const std::vector<double>& zPosMeasure,
                          const std::vector<double>& xResolMeasure,
                          const std::vector<double>& yResolMeasure,
                          double& xOriginLine,
                          double& yOriginLine,
                          double& xAngleLine,
                          double& yAngleLine,
                          double& xChisqLine,
                          double& yChisqLine,
                          std::vector<double>& xResidMeasure,
                          std::vector<double>& yResidMeasure);


static const std::string help_usage = R"(
Usage:
  -help                    help message
  -verbose                 verbose flag
  -hitFile         [PATH]   path to hit file  (input)
  -geometryFile    [PATH]   geometry input file (input)
  -rootFile        [PATH]   path to root file (output)
  -maxEventNumber   [int]    max number of events to be processed
  -maxTrackNumber   [int]    max number of tracks to be processed

examples:
../bin/LineFit_main -hitFile /work/data/TB2006/alpide_200629033515.json  -geo align_313_geo.json -r linefit.root -maxT 1000
)";

int main(int argc, char *argv[]) {
  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"geometryFile", required_argument, NULL, 'g'},
                              {"hitFile", required_argument, NULL, 'o'},
                              {"rootFile", required_argument, NULL, 'u'},
                              {"maxEventNumber", required_argument, NULL, 'm'},
                              {"maxTrackNumber", required_argument, NULL, 'n'},
                             {0, 0, 0, 0}};

  std::string geometryFile_path;
  std::string hitFile_path;
  std::string rootFile_path;
  size_t maxTrackNumber = -1;
  size_t maxEventNumber = -1;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
    case 'g':
      geometryFile_path = optarg;
      break;
    case 'o':
      hitFile_path = optarg;
      break;
    case 'u':
      rootFile_path = optarg;
      break;
    case 'm':
      maxEventNumber = std::stoull(optarg);
      break;
    case 'n':
      maxTrackNumber = std::stoull(optarg);
      break;
      /////generic part below///////////
    case 0: /* getopt_long() set a variable, just keep going */
      break;
    case 1:
      fprintf(stderr, "case 1\n");
      exit(1);
      break;
    case ':':
      fprintf(stderr, "case :\n");
      exit(1);
      break;
    case '?':
      fprintf(stderr, "case ?\n");
      exit(1);
      break;
    default:
      fprintf(stderr, "case default, missing branch in switch-case\n");
      exit(1);
      break;
    }
  }
  if(do_help){
    std::fprintf(stdout, "help message:\n");
    std::fprintf(stdout, "%s\n", help_usage.c_str());
    std::exit(1);
  }

  if (geometryFile_path.empty() ||
      hitFile_path.empty()||
      rootFile_path.empty()
    ) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    std::exit(0);
  }

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "geometryFile: %s\n", geometryFile_path.c_str());
  std::fprintf(stdout, "hitFile:      %s\n", hitFile_path.c_str());
  std::fprintf(stdout, "rootFile:     %s\n", rootFile_path.c_str());
  std::fprintf(stdout, "\n");

  std::printf("--------read geo-----\n");
  std::string str_geo = JsonUtils::readFile(geometryFile_path.c_str());
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", geometryFile_path.c_str() );
    throw;
  }

  std::map<size_t, Eigen::Affine3d> trafo_targets;
  if(jsd_geo.HasMember("geometry")&&
     jsd_geo["geometry"].HasMember("targets") &&
     jsd_geo["geometry"]["targets"].IsArray()){
    const auto &js_duts = jsd_geo["geometry"]["targets"];
    for(const auto& js_dut : js_duts.GetArray()){
      size_t id = js_dut["id"].GetUint();
      double cx = js_dut["center"]["x"].GetDouble();
      double cy = js_dut["center"]["y"].GetDouble();
      double cz = js_dut["center"]["z"].GetDouble();
      double rx = js_dut["rotation"]["x"].GetDouble();
      double ry = js_dut["rotation"]["y"].GetDouble();
      double rz = js_dut["rotation"]["z"].GetDouble();

      Eigen::AngleAxisd rotZ(rz, Eigen::Vector3d::UnitZ());
      Eigen::AngleAxisd rotY(ry, Eigen::Vector3d::UnitY());
      Eigen::AngleAxisd rotX(rx, Eigen::Vector3d::UnitX());
      Eigen::Affine3d trafo(Eigen::Affine3d::Identity());
      trafo.rotate(rotX);
      trafo.rotate(rotY);
      trafo.rotate(rotZ);
      trafo.translation() =  Eigen::Vector3d(cx, cy, cz);
      trafo_targets[id]=trafo;
    }
  }

  if( !jsd_geo.HasMember("geometry") ||
      !jsd_geo["geometry"].HasMember("detectors") ||
      !jsd_geo["geometry"]["detectors"].IsArray()){
    std::fprintf(stderr, "Geometry file <%s>: unknown json structure.\n", geometryFile_path.c_str() );
    JsonUtils::printJsonValue(jsd_geo, true);
    throw;
  }
  const auto &js_geo = jsd_geo["geometry"];
  const auto &js_dets = js_geo["detectors"];
  for(const auto& js_det : js_dets.GetArray()){
    JsonUtils::printJsonValue(js_det, false);
  }

  std::map<size_t, Eigen::Affine3d> trafo_dets;
  for(const auto& js_det : js_dets.GetArray()){
    size_t id = js_det["id"].GetUint();
    double cx = js_det["center"]["x"].GetDouble();
    double cy = js_det["center"]["y"].GetDouble();
    double cz = js_det["center"]["z"].GetDouble();
    double rx = js_det["rotation"]["x"].GetDouble();
    double ry = js_det["rotation"]["y"].GetDouble();
    double rz = js_det["rotation"]["z"].GetDouble();

    Eigen::AngleAxisd rotZ(rz, Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd rotY(ry, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd rotX(rx, Eigen::Vector3d::UnitX());
    Eigen::Affine3d trafoGlobalFromMeas(Eigen::Affine3d::Identity());
    trafoGlobalFromMeas.rotate(rotX);
    trafoGlobalFromMeas.rotate(rotY);
    trafoGlobalFromMeas.rotate(rotZ);
    trafoGlobalFromMeas.translation() =  Eigen::Vector3d(cx, cy, cz);
    trafo_dets[id]=trafoGlobalFromMeas;
  }

  TTree tree("linefit", "linefit");
  size_t n_datapack_select_opt = 20000;
  JsonFileDeserializer jsf(hitFile_path);

  std::vector<size_t> idMeas;
  std::vector<double> xPosMeas;
  std::vector<double> yPosMeas;
  std::vector<double> zPosMeas;
  std::vector<double> xResolMeas;
  std::vector<double> yResolMeas;

  std::vector<double> xResidFit;
  std::vector<double> yResidFit;
  std::vector<double> xPosFit;
  std::vector<double> yPosFit;

  std::vector<size_t>* p_idMeas = &idMeas;
  std::vector<double>* p_xPosMeas = &xPosMeas;
  std::vector<double>* p_yPosMeas = &yPosMeas;
  std::vector<double>* p_zPosMeas = &zPosMeas;
  std::vector<double>* p_xResolMeas = &xResolMeas;
  std::vector<double>* p_yResolMeas = &yResolMeas;
  std::vector<double>* p_xResidFit = &xResidFit;
  std::vector<double>* p_yResidFit = &yResidFit;
  std::vector<double>* p_xPosFit = &xPosFit;
  std::vector<double>* p_yPosFit = &yPosFit;

  auto brId = tree.Branch("id", &p_idMeas);
  auto brXMeas = tree.Branch("xmeas", &p_xPosMeas);
  auto brYMeas = tree.Branch("ymeas", &p_yPosMeas);
  auto brXResid = tree.Branch("xresid", &p_xResidFit);
  auto brYResid = tree.Branch("yresid", &p_yResidFit);
  auto brXFit = tree.Branch("xfit", &p_xPosFit);
  auto brYFit = tree.Branch("yfit", &p_yPosFit);



  std::vector<size_t> idTarget;
  std::vector<double> xMeasTarget;
  std::vector<double> yMeasTarget;
  std::vector<double> xFitMeasTarget;
  std::vector<double> yFitMeasTarget;
  std::vector<double> xResidMeasTarget;
  std::vector<double> yResidMeasTarget;

  std::vector<size_t>* p_idTarget = &idTarget;
  std::vector<double>* p_xMeasTarget = &xMeasTarget;
  std::vector<double>* p_yMeasTarget = &yMeasTarget;
  std::vector<double>* p_xFitMeasTarget = &xFitMeasTarget;
  std::vector<double>* p_yFitMeasTarget = &yFitMeasTarget;
  std::vector<double>* p_xResidMeasTarget = &xResidMeasTarget;
  std::vector<double>* p_yResidMeasTarget = &yResidMeasTarget;

  auto brIdTarget = tree.Branch("idTarget", &p_idTarget);
  auto brXMeasTarget = tree.Branch("xMeasTarget", &p_xMeasTarget);
  auto brYMeasTarget = tree.Branch("yMeasTarget", &p_yMeasTarget);
  auto brXFitMeasTarget = tree.Branch("xFitMeasTarget", &p_xFitMeasTarget);
  auto brYFitMeasTarget = tree.Branch("yFitMeasTarget", &p_yFitMeasTarget);
  auto brXResidMeasTarget = tree.Branch("xResidMeasTarget", &p_xResidMeasTarget);
  auto brYResidMeasTarget = tree.Branch("yResidMeasTarget", &p_yResidMeasTarget);

  size_t nTracks = 0;
  size_t nEvents = 0;
  while (jsf) {
    if (nTracks >= maxTrackNumber || nEvents >= maxEventNumber ) {
      break;
    }
    auto evpack = jsf.getNextJsonDocument();
    if(evpack.IsNull()){
      std::fprintf(stdout, "reach null object, end of file\n");
      break;
    }
    const auto &layers = evpack["layers"];

    bool isAllSingleHit = true;
    for (const auto &layer : layers.GetArray()){
      if(layer["hit"].Size() !=1){
        isAllSingleHit = false;
        break;
      }
    }
    if(!isAllSingleHit){
      std::fprintf(stdout, "isAllSingleHit false, skip envent %i\n", nEvents);
      nEvents++;
      continue;
    }

    idMeas.clear();
    xPosMeas.clear();
    yPosMeas.clear();
    zPosMeas.clear();
    xResolMeas.clear();
    yResolMeas.clear();

    xResidFit.clear();
    yResidFit.clear();
    xPosFit.clear();
    yPosFit.clear();

    for(const auto& [theId, theTrafo]: trafo_dets){
      for (const auto &layer : layers.GetArray()){
        size_t id = layer["ext"].GetUint();
        if(id != theId){
          continue;
        }
        if(layer["hit"].Size() !=1){
          throw;
        }
        const auto &hit = layer["hit"][0];
        double x_hit = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
        double y_hit = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;
        Eigen::Vector3d measPos(x_hit, y_hit, 0);
        idMeas.push_back(id);
        Eigen::Vector3d globalPos =  theTrafo * measPos;
        xPosMeas.push_back(globalPos(0));
        yPosMeas.push_back(globalPos(1));
        zPosMeas.push_back(globalPos(2));
        xResolMeas.push_back(0.02);
        yResolMeas.push_back(0.02);
      }
    }

//
    double xOriLine;
    double yOriLine;
    double xAngleLine;
    double yAngleLine;
    double xChisqLine;
    double yChisqLine;
    FitTrack(idMeas.size(), xPosMeas, yPosMeas, zPosMeas, xResolMeas, yResolMeas,
             xOriLine, yOriLine, xAngleLine, yAngleLine, xChisqLine, yChisqLine,
             xResidFit, yResidFit);

    for(int i=0; i<idMeas.size(); i++){
      xPosFit.push_back(xPosMeas[i]-xResidFit[i]);
      yPosFit.push_back(yPosMeas[i]-yResidFit[i]);
    }

    xMeasTarget.clear();
    yMeasTarget.clear();
    xFitMeasTarget.clear();
    yFitMeasTarget.clear();
    xResidMeasTarget.clear();
    yResidMeasTarget.clear();

    Eigen::Vector3d oriLine(xOriLine, yOriLine, 0.);
    Eigen::Vector3d dirLine(tan(xAngleLine), tan(yAngleLine), 1);
    dirLine.normalize();
    Eigen::ParametrizedLine<double, 3> line(oriLine, dirLine);

    for(const auto& [theId, theTrafo]: trafo_targets){
      for (const auto &layer : layers.GetArray()){
        size_t id = layer["ext"].GetUint();
        if(id != theId){
          continue;
        }
        if(layer["hit"].Size() != 1){
          throw;
        }
        const auto &hit = layer["hit"][0];
        double x_hit = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
        double y_hit = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;
        Eigen::Vector3d measPos(x_hit, y_hit, 0);
        idTarget.push_back(id);
        xMeasTarget.push_back(measPos(0));
        yMeasTarget.push_back(measPos(1));

        Eigen::Vector3d ori = theTrafo*Eigen::Vector3d(0., 0., 0.);
        Eigen::Vector3d dir = theTrafo*Eigen::Vector3d(0., 0., 1.);
        Eigen::Hyperplane<double, 3> targetPlane(dir, ori);
        Eigen::Vector3d fitGlobal = line.intersectionPoint(targetPlane);
        Eigen::Vector3d fitLocal = theTrafo.inverse() * fitGlobal;
        xFitMeasTarget.push_back(fitLocal(0));
        yFitMeasTarget.push_back(fitLocal(1));
        Eigen::Vector3d residMeas = measPos - fitLocal;
        xResidMeasTarget.push_back(residMeas(0));
        yResidMeasTarget.push_back(residMeas(1));
      }
    }

    tree.Fill();
    nEvents++;
    nTracks++;
  }

  TFile tfile(rootFile_path.c_str(),"recreate");
  tree.Write();
  tfile.Close();
  return 0;
}


void FitTrack(unsigned int nMeasures,
              const std::vector<double>& xPosMeasure,
              const std::vector<double>& yPosMeasure,
              const std::vector<double>& zPosMeasure,
              const std::vector<double>& xResolMeasure,
              const std::vector<double>& yResolMeasure,
              double& xOriginLine,
              double& yOriginLine,
              double& xAngleLine,
              double& yAngleLine,
              double& xChisqLine,
              double& yChisqLine,
              std::vector<double>& xResidMeasure,
              std::vector<double>& yResidMeasure){
  unsigned int nPlanesFit = nMeasures;
  const double * xPosFit = xPosMeasure.data();
  const double * yPosFit = yPosMeasure.data();
  const double * zPosFit = zPosMeasure.data();;
  const double * xResFit = xResolMeasure.data();
  const double * yResFit = yResolMeasure.data();

  xResidMeasure.resize(nPlanesFit);
  yResidMeasure.resize(nPlanesFit);
  double *residXFit = xResidMeasure.data();
  double *residYFit = yResidMeasure.data();

  // std::cout<< "------------"<<std::endl;
  // std::cout<< nMeasures<<std::endl;
  // std::cout<< xPosMeasure.size()<<std::endl;
  // std::cout<< yPosMeasure.size()<<std::endl;
  // std::cout<< zPosMeasure.size()<<std::endl;
  // std::cout<< xResolMeasure.size()<<std::endl;
  // std::cout<< yResolMeasure.size()<<std::endl;
  // std::cout<< xResidMeasure.size()<<std::endl;
  // std::cout<< yResidMeasure.size()<<std::endl;
  // std::cout<< "============="<<std::endl;


  float S1[2]   = {0,0};
  float Sx[2]   = {0,0};
  float Xbar[2] = {0,0};

  std::vector<float> Zbar_X_v(nPlanesFit);
  std::vector<float> Zbar_Y_v(nPlanesFit);

  float * Zbar_X = Zbar_X_v.data();
  float * Zbar_Y = Zbar_Y_v.data();

  for (unsigned int counter = 0; counter < nPlanesFit; counter++){
    Zbar_X[counter] = 0.;
    Zbar_Y[counter] = 0.;
  }

  float Sy[2]     = {0,0};
  float Ybar[2]   = {0,0};
  float Sxybar[2] = {0,0};
  float Sxxbar[2] = {0,0};
  float A2[2]     = {0,0};

  // define S1
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    S1[0] = S1[0] + 1/pow(xResFit[counter],2);
    S1[1] = S1[1] + 1/pow(yResFit[counter],2);
  }

  // define Sx
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    Sx[0] = Sx[0] + zPosFit[counter]/pow(xResFit[counter],2);
    Sx[1] = Sx[1] + zPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Xbar
  Xbar[0]=Sx[0]/S1[0];
  Xbar[1]=Sx[1]/S1[1];

  // coordinate transformation !! -> bar
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    Zbar_X[counter] = zPosFit[counter]-Xbar[0];
    Zbar_Y[counter] = zPosFit[counter]-Xbar[1];
  }

  // define Sy
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    Sy[0] = Sy[0] + xPosFit[counter]/pow(xResFit[counter],2);
    Sy[1] = Sy[1] + yPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Ybar
  Ybar[0]=Sy[0]/S1[0];
  Ybar[1]=Sy[1]/S1[1];

  // define Sxybar
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    Sxybar[0] = Sxybar[0] + Zbar_X[counter] * xPosFit[counter]/pow(xResFit[counter],2);
    Sxybar[1] = Sxybar[1] + Zbar_Y[counter] * yPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Sxxbar
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    Sxxbar[0] = Sxxbar[0] + Zbar_X[counter] * Zbar_X[counter]/pow(xResFit[counter],2);
    Sxxbar[1] = Sxxbar[1] + Zbar_Y[counter] * Zbar_Y[counter]/pow(yResFit[counter],2);
  }

  // define A2
  A2[0]=Sxybar[0]/Sxxbar[0];
  A2[1]=Sxybar[1]/Sxxbar[1];

  // Calculate chi sqaured
  // Chi^2 for X and Y coordinate for hits in all planes
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    xChisqLine += pow(  -zPosFit[counter]*A2[0] +  xPosFit[counter] - Ybar[0] + Xbar[0]*A2[0]  ,2) /
      pow(  xResFit[counter]  ,2);

    xChisqLine += pow(-zPosFit[counter]*A2[1]
                      +yPosFit[counter]-Ybar[1]+Xbar[1]*A2[1],2)/pow(yResFit[counter],2);
  }

  for(unsigned int counter = 0; counter < nPlanesFit; counter++ ) {
    residXFit[counter] = xPosFit[counter] - (Ybar[0]-Xbar[0]*A2[0]+zPosFit[counter]*A2[0]); // sign reverse?
    residYFit[counter] = yPosFit[counter] - (Ybar[1]-Xbar[1]*A2[1]+zPosFit[counter]*A2[1]);
  }

  // define angle
  xAngleLine = atan(A2[0]);
  yAngleLine = atan(A2[1]);

  xOriginLine = Ybar[0] - Xbar[0]*A2[0];
  yOriginLine = Ybar[1] - Xbar[1]*A2[1];
}
