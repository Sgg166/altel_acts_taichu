// origin from: Igor Rubinskiy
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <iostream>
#include <fstream>

#include <cstring>
#include <cmath>

#include <Eigen/Geometry>

#include "Mille.h"
#include "TelMille.hh"
#include "exampleUtil.h"

using namespace std;

EUTelMille::EUTelMille(){

}


EUTelMille::~EUTelMille(){
};


/// Create a silicon layer with 2D measurement.
/**
 * Create silicon layer with 2D measurement (u,v) at fixed Z-position.
 * The measurement directions in the XY plane can be orthogonal or non-orthogonal
 * (but must be different).
 *
 * \param [in] aName      name
 * \param [in] layer      layer ID
 * \param [in] xPos       X-position (of center)
 * \param [in] yPos       Y-position (of center)
 * \param [in] zPos       Z-position (of center)
 * \param [in] thickness  thickness / radiation_length
 * \param [in] uAngle     angle of u-direction in XY plane
 * \param [in] uRes       resolution in u-direction
 * \param [in] vAngle     angle of v-direction in XY plane
 * \param [in] vRes       resolution in v-direction
 */
std::unique_ptr<gbl::GblDetectorLayer> EUTelMille::CreateLayerSit(const std::string& aName, unsigned int layer,
                                                                  double xPos, double yPos, double zPos, double thickness,
                                                                  double uAngle, double uRes, double vAngle, double vRes) {
  Eigen::Vector3d aCenter(xPos, yPos, zPos);
  Eigen::Vector2d aResolution(uRes, vRes);
  Eigen::Vector2d aPrecision(1. / (uRes * uRes), 1. / (vRes * vRes));
  Eigen::Matrix3d measTrafo;

  const double cosU = cos(uAngle / 180. * M_PI);
  const double sinU = sin(uAngle / 180. * M_PI);

  const double cosV = cos(vAngle / 180. * M_PI);
  const double sinV = sin(vAngle / 180. * M_PI);

  measTrafo << cosU, sinU, 0.,   cosV, sinV, 0.,     0., 0., 1.; // U,V,N
  Eigen::Matrix3d alignTrafo;
  alignTrafo << 1., 0., 0.,      0., 1., 0.,         0., 0., 1.; // X,Y,Z
  return std::make_unique<gbl::GblDetectorLayer>(aName, layer, 2, thickness, aCenter, aResolution,
                                                 aPrecision, measTrafo, alignTrafo);
}




/*! Performs analytic straight line fit.
 *
 * Determines parameters of a straight line passing through the
 * measurement points.
 *
 * @see http://www.desy.de/~blobel/eBuch.pdf page 162
 *
 */
void EUTelMille::FitTrack(unsigned int nMeasures,
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
                          std::vector<double>& yResidMeasure
  ) {
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

void EUTelMille::setResolution(double resolX, double resolY){
  for(auto &[id, n]: m_indexDet){
    m_xResolution[id] = resolX;
    m_yResolution[id] = resolY;
  }
}

void EUTelMille::setResolution(size_t id, double resolX, double resolY){
  m_xResolution[id] = resolX;
  m_yResolution[id] = resolY;
}


void EUTelMille::setGeometry(const JsonValue& js) {
  if(!js.HasMember("geometry")){
    std::fprintf(stderr, "unable to find \"geomerty\" key for detector geomerty from JS\n");
    throw;
  }

  const auto &js_geo = js["geometry"];
  const auto &js_dets = js_geo["detectors"];
  std::map<double, size_t> zmap_sort;
  for(const auto& js_det: js_dets.GetArray()){
    size_t id = js_det["id"].GetUint();
    double cz = js_det["center"]["z"].GetDouble();
    zmap_sort[cz]=id;
  }

  m_nPlanes=zmap_sort.size();
  size_t n=0;
  for(auto [the_cz, the_id]: zmap_sort){
    for(const auto& js_det: js_dets.GetArray()){
      size_t id = js_det["id"].GetUint();
      if(id == the_id){
        double cx = js_det["center"]["x"].GetDouble();
        double cy = js_det["center"]["y"].GetDouble();
        double cz = js_det["center"]["z"].GetDouble();
        double rx = js_det["rotation"]["x"].GetDouble();
        double ry = js_det["rotation"]["y"].GetDouble();
        double rz = js_det["rotation"]["z"].GetDouble();
        m_xPosDet[id]=cx;
        m_yPosDet[id]=cy;
        m_zPosDet[id]=cz;

        m_alphaPosDet[id]=rx;
        m_betaPosDet[id]=ry;
        m_gammaPosDet[id]=rz;
        m_indexDet[id] = n;
        m_dets[id]=CreateLayerSit("PIX"+std::to_string(id), id,
                                  cx, cy, cz, 0.001,
                                  rz, 0.02, rz+90., 0.02);

        n++;
      }
    }
  }
}

void EUTelMille::startMilleBinary(const std::string& path){
  m_mille.reset(new Mille(path.c_str()));
  m_binPath=path;
}


void EUTelMille::endMilleBinary(){
  m_mille.reset();
}

void EUTelMille::fillTrackXYRz(const JsonValue& js) {
  double stdevFact = 0.09;

  std::vector<double> xPosHit(m_nPlanes,0);
  std::vector<double> yPosHit(m_nPlanes,0);
  std::vector<double> zPosHit(m_nPlanes,0);

  std::vector<double> xResolHit(m_nPlanes, 0);
  std::vector<double> yResolHit(m_nPlanes, 0);

  std::vector<size_t> idHit(m_nPlanes,0);

  if(js.Size()!= m_nPlanes){
    std::fprintf(stderr, "hits number[%i] is less than detector number[%i] \n", js.Size(), m_nPlanes);
    throw;
  }
  for(const auto& js_hit : js.GetArray()){
    size_t id =js_hit["id"].GetUint();
    size_t detN = m_indexDet.at(id);
    double xPosDet = m_xPosDet.at(id);
    double yPosDet = m_yPosDet.at(id);
    double zPosDet = m_zPosDet.at(id);
    double xRotDet = m_alphaPosDet.at(id);
    double yRotDet = m_betaPosDet.at(id);
    double zRotDet = m_gammaPosDet.at(id);

    Eigen::AngleAxisd rotZ(zRotDet, Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd rotY(yRotDet, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd rotX(xRotDet, Eigen::Vector3d::UnitX());
    Eigen::Affine3d trafoGlobalFromMeas(Eigen::Affine3d::Identity());
    trafoGlobalFromMeas.rotate(rotX);
    trafoGlobalFromMeas.rotate(rotY);
    trafoGlobalFromMeas.rotate(rotZ);
    trafoGlobalFromMeas.translation() =  Eigen::Vector3d(xPosDet, yPosDet, zPosDet);

    Eigen::Vector3d measPos(js_hit["x"].GetDouble(), js_hit["y"].GetDouble(), 0);
    Eigen::Vector3d globalPos =  trafoGlobalFromMeas * measPos;

    xPosHit[detN]=globalPos[0];
    yPosHit[detN]=globalPos[1];
    zPosHit[detN]=globalPos[2];

    //TODO: from meas UV to XY
    xResolHit[detN]=m_xResolution.at(id);
    yResolHit[detN]=m_yResolution.at(id);

    idHit[detN]=id;
  }


  double xOriginTrack = 0;
  double yOriginTrack = 0;
  double xAngleTrack = 0;
  double yAngleTrack = 0;
  double xChisqTrack = 0;
  double yChisqTrack = 0;

  std::vector<double> xResidHit(m_nPlanes,0);
  std::vector<double> yResidHit(m_nPlanes,0);

  // Calculate residuals
  FitTrack(m_nPlanes,
           xPosHit,
           yPosHit,
           zPosHit,
           xResolHit,
           yResolHit,

           xOriginTrack,
           yOriginTrack,
           xAngleTrack,
           yAngleTrack,

           xChisqTrack,
           yChisqTrack,

           xResidHit,
           yResidHit
    );

  const int nLC = 4; // number of local parameters, x, y, xa, ya
  const int nGL = m_nPlanes * 3; // number of global parameters, x, y, rz

  std::vector<float> derLC(nLC, 0);
  std::vector<float> derGL(nGL, 0);

  std::vector<int> label(nGL, 0);

  for (unsigned int n = 0; n < m_nPlanes; n++) {
    unsigned int id=-1;
    for(auto [the_id, the_n]: m_indexDet){
      if(the_n == n){
        id = the_id;
        break;
      }
    }
    if(id == -1){
      std::cout<< "not find id by n number"<<std::endl;
      throw;
    }
    auto it_label_layer_begin = std::begin(label);
    std::advance(it_label_layer_begin, 3*n);
    auto it_label_layer_end   = std::begin(label);
    std::advance(it_label_layer_end, 3*(n+1));
    std::iota(it_label_layer_begin, it_label_layer_end, id*10+1);
  }

  // std::iota(std::begin(label), std::end(label), 1);

  // loop over all planes
  for (unsigned int n = 0; n < m_nPlanes; n++) {
    unsigned int id=-1;
    for(auto [the_id, the_n]: m_indexDet){
      if(the_n == n){
        id = the_id;
        break;
      }
    }
    if(id == -1){
      std::cout<< "not find id by n number"<<std::endl;
      throw;
    }
    Eigen::Vector3d posPredit( xPosHit[n]-xResidHit[n], yPosHit[n]-yResidHit[n], zPosHit[n]);
    Eigen::Vector3d dirLine( tan(xAngleTrack), tan(yAngleTrack), 1);

    auto& det = m_dets.at(id);
    Eigen::Matrix<double, 2, 6> fullDerGL = det->getRigidBodyDerLocal_mod(posPredit, dirLine); // global
// or global cordination?


    float residual;
    float sigma;

    // derGL[((n * 3) + 0)] = -1;
    // derGL[((n * 3) + 2)] = yPosHit[n];
    derGL[((n * 3) + 0)] = fullDerGL(0,0); // 1
    derGL[((n * 3) + 1)] = fullDerGL(0,1); // 0
    derGL[((n * 3) + 2)] = fullDerGL(0,5); // y

    derLC[0] = 1;
    derLC[2] = zPosHit[n];
    residual = xResidHit[n];
    sigma    = xResolHit[n];
    m_mille->mille(nLC,derLC.data(),nGL,derGL.data(),label.data(),residual,sigma);

    derGL[((n * 3) + 0)] = 0;
    derGL[((n * 3) + 1)] = 0;
    derGL[((n * 3) + 2)] = 0;
    derLC[0] = 0;
    derLC[2] = 0;

    // derGL[((n * 3) + 1)] = -1;
    // derGL[((n * 3) + 2)] = -1 * xPosHit[n];
    derGL[((n * 3) + 0)] = fullDerGL(1,0); // 0
    derGL[((n * 3) + 1)] = fullDerGL(1,1); // 1
    derGL[((n * 3) + 2)] = fullDerGL(1,5); // -x


    derLC[1] = 1;
    derLC[3] = zPosHit[n];
    residual = yResidHit[n];
    sigma    = yResolHit[n];
    m_mille->mille(nLC,derLC.data(),nGL,derGL.data(),label.data(),residual,sigma);

    derGL[((n * 3) + 0)] = 0;
    derGL[((n * 3) + 1)] = 0;
    derGL[((n * 3) + 2)] = 0;
    derLC[1] = 0;
    derLC[3] = 0;

  } // end loop over all planes

  m_mille->end();
}


void EUTelMille::createPedeStreeringModeXYRz(const std::string& path){
  ofstream steerFile;
  steerFile.open(path.c_str());

  steerFile << "Cfiles" << endl;
  steerFile << m_binPath << endl;
  steerFile << endl;

  steerFile << "Parameter" << endl;

  // first (minimium z) plane are all fixed; last plane is center fixed, but rotatable
  for (unsigned int n = 0; n < m_nPlanes; n++) {
    for(auto &[id, detN]: m_indexDet ){
      if(detN == n){
        if(n == 0){
          steerFile << (id*10 + 1) << " " << "0.0" << " -1.0" << endl;
          steerFile << (id*10 + 2) << " " << "0.0" << " -1.0" << endl;
          steerFile << (id*10 + 3) << " " << "0.0" << " -1.0" << endl;
        }
        else if(n == m_nPlanes-1){
          steerFile << (id*10 + 1) << " " << "0.0" << " -1.0" << endl;
          steerFile << (id*10 + 2) << " " << "0.0" << " -1.0" << endl;
          steerFile << (id*10 + 3) << " " << "0.0" << " 0.0" << endl;
        }
        else{
          steerFile << (id*10 + 1) << " " << "0.0" << " 0.0" << endl;
          steerFile << (id*10 + 2) << " " << "0.0" << " 0.0" << endl;
          steerFile << (id*10 + 3) << " " << "0.0" << " 0.0" << endl;
        }
        break;
      }
    }
  }

  steerFile << endl;
  steerFile << endl;
  steerFile << "method inversion 10 0.0001" << endl;
  steerFile << endl;
  steerFile << "histprint" << endl;
  steerFile << endl;
  steerFile << "end" << endl;
  steerFile.close();
}
