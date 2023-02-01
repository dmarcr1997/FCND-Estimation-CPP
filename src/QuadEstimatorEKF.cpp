#include "Common.h"
#include "QuadEstimatorEKF.h"
#include "Utility/SimpleConfig.h"
#include "Utility/StringUtils.h"
#include "Math/Quaternion.h"

using namespace SLR;

const int QuadEstimatorEKF::QUAD_EKF_NUM_STATES;

QuadEstimatorEKF::QuadEstimatorEKF(string config, string name)
  : BaseQuadEstimator(config),
  Q(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES),
  R_GPS(6, 6),
  R_Mag(1, 1),
  ekfState(QUAD_EKF_NUM_STATES),
  ekfCov(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES),
  trueError(QUAD_EKF_NUM_STATES)
{
  _name = name;
  Init();
}

QuadEstimatorEKF::~QuadEstimatorEKF()
{

}

void QuadEstimatorEKF::Init()
{
  ParamsHandle paramSys = SimpleConfig::GetInstance();

  paramSys->GetFloatVector(_config + ".InitState", ekfState);

  VectorXf initStdDevs(QUAD_EKF_NUM_STATES);
  paramSys->GetFloatVector(_config + ".InitStdDevs", initStdDevs);
  ekfCov.setIdentity();
  for (int i = 0; i < QUAD_EKF_NUM_STATES; i++)
  {
    ekfCov(i, i) = initStdDevs(i) * initStdDevs(i);
  }

  // complementary filter params
  attitudeTau = paramSys->Get(_config + ".AttitudeTau", .1f);
  dtIMU = paramSys->Get(_config + ".dtIMU", .002f);

  pitchEst = 0;
  rollEst = 0;
  
  // GPS measurement model covariance
  R_GPS.setZero();
  R_GPS(0, 0) = R_GPS(1, 1) = powf(paramSys->Get(_config + ".GPSPosXYStd", 0), 2);
  R_GPS(2, 2) = powf(paramSys->Get(_config + ".GPSPosZStd", 0), 2);
  R_GPS(3, 3) = R_GPS(4, 4) = powf(paramSys->Get(_config + ".GPSVelXYStd", 0), 2);
  R_GPS(5, 5) = powf(paramSys->Get(_config + ".GPSVelZStd", 0), 2);

  // magnetometer measurement model covariance
  R_Mag.setZero();
  R_Mag(0, 0) = powf(paramSys->Get(_config + ".MagYawStd", 0), 2);

  // load the transition model covariance
  Q.setZero();
  Q(0, 0) = Q(1, 1) = powf(paramSys->Get(_config + ".QPosXYStd", 0), 2);
  Q(2, 2) = powf(paramSys->Get(_config + ".QPosZStd", 0), 2);
  Q(3, 3) = Q(4, 4) = powf(paramSys->Get(_config + ".QVelXYStd", 0), 2);
  Q(5, 5) = powf(paramSys->Get(_config + ".QVelZStd", 0), 2);
  Q(6, 6) = powf(paramSys->Get(_config + ".QYawStd", 0), 2);
  Q *= dtIMU;

  rollErr = pitchErr = maxEuler = 0;
  posErrorMag = velErrorMag = 0;
}

void QuadEstimatorEKF::UpdateFromIMU(V3F accel, V3F gyro)
{
  float phi = rollEst;
  float theta = pitchEst;

  Mat3x3F rotationMatrix = Mat3x3F::Zeros();
  rotationMatrix(0, 0) = 1;
  rotationMatrix(0, 1) = sin(phi) * tan(theta);
  rotationMatrix(0, 2) = cos(phi) * tan(theta);
  rotationMatrix(1, 1) = cos(phi);
  rotationMatrix(1, 2) = -sin(phi);
  rotationMatrix(2, 1) = sin(phi) / cos(theta);
  rotationMatrix(2, 2) = cos(phi) / cos(theta);

  V3F angleDot = rotationMatrix * gyro;

  float predictedPitch = pitchEst + dtIMU * angleDot.y;
  float predictedRoll = rollEst + dtIMU * angleDot.x;
  ekfState(6) = ekfState(6) + dtIMU * angleDot.z;	// yaw

  // normalize yaw to -pi .. pi
  if (ekfState(6) > F_PI) ekfState(6) -= 2.f*F_PI;
  if (ekfState(6) < -F_PI) ekfState(6) += 2.f*F_PI;

  // CALCULATE UPDATE
  accelRoll = atan2f(accel.y, accel.z);
  accelPitch = atan2f(-accel.x, 9.81f);

  // FUSE INTEGRATION AND UPDATE
  rollEst = attitudeTau / (attitudeTau + dtIMU) * (predictedRoll)+dtIMU / (attitudeTau + dtIMU) * accelRoll;
  pitchEst = attitudeTau / (attitudeTau + dtIMU) * (predictedPitch)+dtIMU / (attitudeTau + dtIMU) * accelPitch;

  lastGyro = gyro;
}

void QuadEstimatorEKF::UpdateTrueError(V3F truePos, V3F trueVel, Quaternion<float> trueAtt)
{
  VectorXf trueState(QUAD_EKF_NUM_STATES);
  trueState(0) = truePos.x;
  trueState(1) = truePos.y;
  trueState(2) = truePos.z;
  trueState(3) = trueVel.x;
  trueState(4) = trueVel.y;
  trueState(5) = trueVel.z;
  trueState(6) = trueAtt.Yaw();

  trueError = ekfState - trueState;
  if (trueError(6) > F_PI) trueError(6) -= 2.f*F_PI;
  if (trueError(6) < -F_PI) trueError(6) += 2.f*F_PI;

  pitchErr = pitchEst - trueAtt.Pitch();
  rollErr = rollEst - trueAtt.Roll();
  maxEuler = MAX(fabs(pitchErr), MAX(fabs(rollErr), fabs(trueError(6))));

  posErrorMag = truePos.dist(V3F(ekfState(0), ekfState(1), ekfState(2)));
  velErrorMag = trueVel.dist(V3F(ekfState(3), ekfState(4), ekfState(5)));
}

VectorXf QuadEstimatorEKF::PredictState(VectorXf curState, float dt, V3F accel, V3F gyro)
{
  assert(curState.size() == QUAD_EKF_NUM_STATES);
  VectorXf predictedState = curState;
  // Predict the current state forward by time dt using current accelerations and body rates as input
  
  Quaternion<float> attitude = Quaternion<float>::FromEuler123_RPY(rollEst, pitchEst, curState(6));

  predictedState(0) = curState(0) + dt * curState(3);
  predictedState(1) = curState(1) + dt * curState(4);
  predictedState(2) = curState(2) + dt * curState(5);

  V3F W = attitude.Rotate_BtoI(accel);

  predictedState(3) = curState(3) + dt * W.x;
  predictedState(4) = curState(4) + dt * W.y;
  predictedState(5) = curState(5) + dt * W.z - dt * CONST_GRAVITY;

  return predictedState;
}

MatrixXf QuadEstimatorEKF::GetRbgPrime(float roll, float pitch, float yaw)
{
  // first, figure out the Rbg_prime
  MatrixXf RbgPrime(3, 3);
  RbgPrime.setZero();

  // Return the partial derivative of the Rbg rotation matrix with respect to yaw. We call this RbgPrime.

  RbgPrime(0, 0) = -cos(pitch) * sin(yaw);
  RbgPrime(0, 1) = -sin(roll) * sin(pitch) * sin(yaw) - cos(roll) * cos(yaw);
  RbgPrime(0, 2) = -cos(roll) * sin(pitch) * sin(yaw) + sin(roll) * cos(yaw);
  RbgPrime(1, 0) = cos(pitch) * cos(yaw);
  RbgPrime(1, 1) = sin(roll) * sin(pitch) * cos(yaw) - cos(roll) * sin(yaw);
  RbgPrime(1, 2) = cos(roll) * sin(pitch) * cos(yaw) + sin(roll) * sin(yaw);
  RbgPrime(2, 0) = 0;
  RbgPrime(2, 1) = 0;
  RbgPrime(2, 2) = 0;

  return RbgPrime;
}

void QuadEstimatorEKF::Predict(float dt, V3F accel, V3F gyro)
{
  // predict the state forward
  VectorXf newState = PredictState(ekfState, dt, accel, gyro);

  // Predict the current covariance forward by dt using the current accelerations and body rates as input.

  // we'll want the partial derivative of the Rbg matrix
  MatrixXf RbgPrime = GetRbgPrime(rollEst, pitchEst, ekfState(6));

  // we've created an empty Jacobian for you, currently simply set to identity
  MatrixXf gPrime(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES);
  gPrime.setIdentity();

  gPrime(0, 3) = dt;
  gPrime(1, 4) = dt;
  gPrime(2, 5) = dt;

  gPrime(3, 6) = (RbgPrime(0) * accel).sum() * dt;
  gPrime(4, 6) = (RbgPrime(1) * accel).sum() * dt;
  gPrime(5, 6) = (RbgPrime(2) * accel).sum() * dt;

  ekfCov = gPrime * ekfCov * gPrime.transpose() + Q;

  ekfState = newState;
}

void QuadEstimatorEKF::UpdateFromGPS(V3F pos, V3F vel)
{
  VectorXf z(6), zFromX(6);
  z(0) = pos.x;
  z(1) = pos.y;
  z(2) = pos.z;
  z(3) = vel.x;
  z(4) = vel.y;
  z(5) = vel.z;

  MatrixXf hPrime(6, QUAD_EKF_NUM_STATES);
  hPrime.setZero();

  zFromX(0) = ekfState(0);
  zFromX(1) = ekfState(1);
  zFromX(2) = ekfState(2);
  zFromX(3) = ekfState(3);
  zFromX(4) = ekfState(4);
  zFromX(5) = ekfState(5);

  for (int i = 0; i < 6; i++) {
      hPrime(i, i) = 1;
  }

  Update(z, hPrime, R_GPS, zFromX);
}

void QuadEstimatorEKF::UpdateFromMag(float magYaw)
{
  VectorXf z(1), zFromX(1);
  z(0) = magYaw;

  MatrixXf hPrime(1, QUAD_EKF_NUM_STATES);
  hPrime.setZero();

  
  zFromX(0) = ekfState(6);
  float difference = magYaw - ekfState(6);
  if (difference > F_PI) {
       zFromX(0) += 2 * F_PI;
  }
  else if (difference < -F_PI) {
       zFromX(0) -= 2 * F_PI;
  }

  hPrime(0, 6) = 1;

  Update(z, hPrime, R_Mag, zFromX);
}

// Execute an EKF update step
// z: measurement
// H: Jacobian of observation function evaluated at the current estimated state
// R: observation error model covariance 
// zFromX: measurement prediction based on current state
void QuadEstimatorEKF::Update(VectorXf& z, MatrixXf& H, MatrixXf& R, VectorXf& zFromX)
{
  assert(z.size() == H.rows());
  assert(QUAD_EKF_NUM_STATES == H.cols());
  assert(z.size() == R.rows());
  assert(z.size() == R.cols());
  assert(z.size() == zFromX.size());

  MatrixXf toInvert(z.size(), z.size());
  toInvert = H*ekfCov*H.transpose() + R;
  MatrixXf K = ekfCov * H.transpose() * toInvert.inverse();

  ekfState = ekfState + K*(z - zFromX);

  MatrixXf eye(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES);
  eye.setIdentity();

  ekfCov = (eye - K*H)*ekfCov;
}

// Calculate the condition number of the EKF ovariance matrix (useful for numerical diagnostics)
// The condition number provides a measure of how similar the magnitudes of the error metric beliefs 
// about the different states are. If the magnitudes are very far apart, numerical issues will start to come up.
float QuadEstimatorEKF::CovConditionNumber() const
{
  MatrixXf m(7, 7);
  for (int i = 0; i < 7; i++)
  {
    for (int j = 0; j < 7; j++)
    {
      m(i, j) = ekfCov(i, j);
    }
  }

  Eigen::JacobiSVD<MatrixXf> svd(m);
  float cond = svd.singularValues()(0)
    / svd.singularValues()(svd.singularValues().size() - 1);
  return cond;
}

// Access functions for graphing variables
bool QuadEstimatorEKF::GetData(const string& name, float& ret) const
{
  if (name.find_first_of(".") == string::npos) return false;
  string leftPart = LeftOf(name, '.');
  string rightPart = RightOf(name, '.');

  if (ToUpper(leftPart) == ToUpper(_name))
  {
#define GETTER_HELPER(A,B) if (SLR::ToUpper(rightPart) == SLR::ToUpper(A)){ ret=(B); return true; }
    GETTER_HELPER("Est.roll", rollEst);
    GETTER_HELPER("Est.pitch", pitchEst);

    GETTER_HELPER("Est.x", ekfState(0));
    GETTER_HELPER("Est.y", ekfState(1));
    GETTER_HELPER("Est.z", ekfState(2));
    GETTER_HELPER("Est.vx", ekfState(3));
    GETTER_HELPER("Est.vy", ekfState(4));
    GETTER_HELPER("Est.vz", ekfState(5));
    GETTER_HELPER("Est.yaw", ekfState(6));

    GETTER_HELPER("Est.S.x", sqrtf(ekfCov(0, 0)));
    GETTER_HELPER("Est.S.y", sqrtf(ekfCov(1, 1)));
    GETTER_HELPER("Est.S.z", sqrtf(ekfCov(2, 2)));
    GETTER_HELPER("Est.S.vx", sqrtf(ekfCov(3, 3)));
    GETTER_HELPER("Est.S.vy", sqrtf(ekfCov(4, 4)));
    GETTER_HELPER("Est.S.vz", sqrtf(ekfCov(5, 5)));
    GETTER_HELPER("Est.S.yaw", sqrtf(ekfCov(6, 6)));

    // diagnostic variables
    GETTER_HELPER("Est.D.AccelPitch", accelPitch);
    GETTER_HELPER("Est.D.AccelRoll", accelRoll);

    GETTER_HELPER("Est.D.ax_g", accelG[0]);
    GETTER_HELPER("Est.D.ay_g", accelG[1]);
    GETTER_HELPER("Est.D.az_g", accelG[2]);

    GETTER_HELPER("Est.E.x", trueError(0));
    GETTER_HELPER("Est.E.y", trueError(1));
    GETTER_HELPER("Est.E.z", trueError(2));
    GETTER_HELPER("Est.E.vx", trueError(3));
    GETTER_HELPER("Est.E.vy", trueError(4));
    GETTER_HELPER("Est.E.vz", trueError(5));
    GETTER_HELPER("Est.E.yaw", trueError(6));
    GETTER_HELPER("Est.E.pitch", pitchErr);
    GETTER_HELPER("Est.E.roll", rollErr);
    GETTER_HELPER("Est.E.MaxEuler", maxEuler);

    GETTER_HELPER("Est.E.pos", posErrorMag);
    GETTER_HELPER("Est.E.vel", velErrorMag);

    GETTER_HELPER("Est.D.covCond", CovConditionNumber());
#undef GETTER_HELPER
  }
  return false;
};

vector<string> QuadEstimatorEKF::GetFields() const
{
  vector<string> ret = BaseQuadEstimator::GetFields();
  ret.push_back(_name + ".Est.roll");
  ret.push_back(_name + ".Est.pitch");

  ret.push_back(_name + ".Est.x");
  ret.push_back(_name + ".Est.y");
  ret.push_back(_name + ".Est.z");
  ret.push_back(_name + ".Est.vx");
  ret.push_back(_name + ".Est.vy");
  ret.push_back(_name + ".Est.vz");
  ret.push_back(_name + ".Est.yaw");

  ret.push_back(_name + ".Est.S.x");
  ret.push_back(_name + ".Est.S.y");
  ret.push_back(_name + ".Est.S.z");
  ret.push_back(_name + ".Est.S.vx");
  ret.push_back(_name + ".Est.S.vy");
  ret.push_back(_name + ".Est.S.vz");
  ret.push_back(_name + ".Est.S.yaw");

  ret.push_back(_name + ".Est.E.x");
  ret.push_back(_name + ".Est.E.y");
  ret.push_back(_name + ".Est.E.z");
  ret.push_back(_name + ".Est.E.vx");
  ret.push_back(_name + ".Est.E.vy");
  ret.push_back(_name + ".Est.E.vz");
  ret.push_back(_name + ".Est.E.yaw");
  ret.push_back(_name + ".Est.E.pitch");
  ret.push_back(_name + ".Est.E.roll");

  ret.push_back(_name + ".Est.E.pos");
  ret.push_back(_name + ".Est.E.vel");

  ret.push_back(_name + ".Est.E.maxEuler");

  ret.push_back(_name + ".Est.D.covCond");

  // diagnostic variables
  ret.push_back(_name + ".Est.D.AccelPitch");
  ret.push_back(_name + ".Est.D.AccelRoll");
  ret.push_back(_name + ".Est.D.ax_g");
  ret.push_back(_name + ".Est.D.ay_g");
  ret.push_back(_name + ".Est.D.az_g");
  return ret;
};
