#include "car.h"

namespace {
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 5.991464547107981;

/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3617894249277739926) {
   out_3617894249277739926[0] = delta_x[0] + nom_x[0];
   out_3617894249277739926[1] = delta_x[1] + nom_x[1];
   out_3617894249277739926[2] = delta_x[2] + nom_x[2];
   out_3617894249277739926[3] = delta_x[3] + nom_x[3];
   out_3617894249277739926[4] = delta_x[4] + nom_x[4];
   out_3617894249277739926[5] = delta_x[5] + nom_x[5];
   out_3617894249277739926[6] = delta_x[6] + nom_x[6];
   out_3617894249277739926[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3153453723933201675) {
   out_3153453723933201675[0] = -nom_x[0] + true_x[0];
   out_3153453723933201675[1] = -nom_x[1] + true_x[1];
   out_3153453723933201675[2] = -nom_x[2] + true_x[2];
   out_3153453723933201675[3] = -nom_x[3] + true_x[3];
   out_3153453723933201675[4] = -nom_x[4] + true_x[4];
   out_3153453723933201675[5] = -nom_x[5] + true_x[5];
   out_3153453723933201675[6] = -nom_x[6] + true_x[6];
   out_3153453723933201675[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_3801013484012494542) {
   out_3801013484012494542[0] = 1.0;
   out_3801013484012494542[1] = 0.0;
   out_3801013484012494542[2] = 0.0;
   out_3801013484012494542[3] = 0.0;
   out_3801013484012494542[4] = 0.0;
   out_3801013484012494542[5] = 0.0;
   out_3801013484012494542[6] = 0.0;
   out_3801013484012494542[7] = 0.0;
   out_3801013484012494542[8] = 0.0;
   out_3801013484012494542[9] = 1.0;
   out_3801013484012494542[10] = 0.0;
   out_3801013484012494542[11] = 0.0;
   out_3801013484012494542[12] = 0.0;
   out_3801013484012494542[13] = 0.0;
   out_3801013484012494542[14] = 0.0;
   out_3801013484012494542[15] = 0.0;
   out_3801013484012494542[16] = 0.0;
   out_3801013484012494542[17] = 0.0;
   out_3801013484012494542[18] = 1.0;
   out_3801013484012494542[19] = 0.0;
   out_3801013484012494542[20] = 0.0;
   out_3801013484012494542[21] = 0.0;
   out_3801013484012494542[22] = 0.0;
   out_3801013484012494542[23] = 0.0;
   out_3801013484012494542[24] = 0.0;
   out_3801013484012494542[25] = 0.0;
   out_3801013484012494542[26] = 0.0;
   out_3801013484012494542[27] = 1.0;
   out_3801013484012494542[28] = 0.0;
   out_3801013484012494542[29] = 0.0;
   out_3801013484012494542[30] = 0.0;
   out_3801013484012494542[31] = 0.0;
   out_3801013484012494542[32] = 0.0;
   out_3801013484012494542[33] = 0.0;
   out_3801013484012494542[34] = 0.0;
   out_3801013484012494542[35] = 0.0;
   out_3801013484012494542[36] = 1.0;
   out_3801013484012494542[37] = 0.0;
   out_3801013484012494542[38] = 0.0;
   out_3801013484012494542[39] = 0.0;
   out_3801013484012494542[40] = 0.0;
   out_3801013484012494542[41] = 0.0;
   out_3801013484012494542[42] = 0.0;
   out_3801013484012494542[43] = 0.0;
   out_3801013484012494542[44] = 0.0;
   out_3801013484012494542[45] = 1.0;
   out_3801013484012494542[46] = 0.0;
   out_3801013484012494542[47] = 0.0;
   out_3801013484012494542[48] = 0.0;
   out_3801013484012494542[49] = 0.0;
   out_3801013484012494542[50] = 0.0;
   out_3801013484012494542[51] = 0.0;
   out_3801013484012494542[52] = 0.0;
   out_3801013484012494542[53] = 0.0;
   out_3801013484012494542[54] = 1.0;
   out_3801013484012494542[55] = 0.0;
   out_3801013484012494542[56] = 0.0;
   out_3801013484012494542[57] = 0.0;
   out_3801013484012494542[58] = 0.0;
   out_3801013484012494542[59] = 0.0;
   out_3801013484012494542[60] = 0.0;
   out_3801013484012494542[61] = 0.0;
   out_3801013484012494542[62] = 0.0;
   out_3801013484012494542[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_7557349604050774031) {
   out_7557349604050774031[0] = state[0];
   out_7557349604050774031[1] = state[1];
   out_7557349604050774031[2] = state[2];
   out_7557349604050774031[3] = state[3];
   out_7557349604050774031[4] = state[4];
   out_7557349604050774031[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_7557349604050774031[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_7557349604050774031[7] = state[7];
}
void F_fun(double *state, double dt, double *out_2200360112784492485) {
   out_2200360112784492485[0] = 1;
   out_2200360112784492485[1] = 0;
   out_2200360112784492485[2] = 0;
   out_2200360112784492485[3] = 0;
   out_2200360112784492485[4] = 0;
   out_2200360112784492485[5] = 0;
   out_2200360112784492485[6] = 0;
   out_2200360112784492485[7] = 0;
   out_2200360112784492485[8] = 0;
   out_2200360112784492485[9] = 1;
   out_2200360112784492485[10] = 0;
   out_2200360112784492485[11] = 0;
   out_2200360112784492485[12] = 0;
   out_2200360112784492485[13] = 0;
   out_2200360112784492485[14] = 0;
   out_2200360112784492485[15] = 0;
   out_2200360112784492485[16] = 0;
   out_2200360112784492485[17] = 0;
   out_2200360112784492485[18] = 1;
   out_2200360112784492485[19] = 0;
   out_2200360112784492485[20] = 0;
   out_2200360112784492485[21] = 0;
   out_2200360112784492485[22] = 0;
   out_2200360112784492485[23] = 0;
   out_2200360112784492485[24] = 0;
   out_2200360112784492485[25] = 0;
   out_2200360112784492485[26] = 0;
   out_2200360112784492485[27] = 1;
   out_2200360112784492485[28] = 0;
   out_2200360112784492485[29] = 0;
   out_2200360112784492485[30] = 0;
   out_2200360112784492485[31] = 0;
   out_2200360112784492485[32] = 0;
   out_2200360112784492485[33] = 0;
   out_2200360112784492485[34] = 0;
   out_2200360112784492485[35] = 0;
   out_2200360112784492485[36] = 1;
   out_2200360112784492485[37] = 0;
   out_2200360112784492485[38] = 0;
   out_2200360112784492485[39] = 0;
   out_2200360112784492485[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2200360112784492485[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2200360112784492485[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2200360112784492485[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2200360112784492485[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2200360112784492485[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2200360112784492485[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2200360112784492485[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2200360112784492485[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2200360112784492485[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2200360112784492485[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2200360112784492485[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2200360112784492485[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2200360112784492485[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2200360112784492485[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2200360112784492485[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2200360112784492485[56] = 0;
   out_2200360112784492485[57] = 0;
   out_2200360112784492485[58] = 0;
   out_2200360112784492485[59] = 0;
   out_2200360112784492485[60] = 0;
   out_2200360112784492485[61] = 0;
   out_2200360112784492485[62] = 0;
   out_2200360112784492485[63] = 1;
}
void h_25(double *state, double *unused, double *out_3119584906358678393) {
   out_3119584906358678393[0] = state[6];
}
void H_25(double *state, double *unused, double *out_185734102743469453) {
   out_185734102743469453[0] = 0;
   out_185734102743469453[1] = 0;
   out_185734102743469453[2] = 0;
   out_185734102743469453[3] = 0;
   out_185734102743469453[4] = 0;
   out_185734102743469453[5] = 0;
   out_185734102743469453[6] = 1;
   out_185734102743469453[7] = 0;
}
void h_24(double *state, double *unused, double *out_5264692651938961980) {
   out_5264692651938961980[0] = state[4];
   out_5264692651938961980[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4282452030868613117) {
   out_4282452030868613117[0] = 0;
   out_4282452030868613117[1] = 0;
   out_4282452030868613117[2] = 0;
   out_4282452030868613117[3] = 0;
   out_4282452030868613117[4] = 1;
   out_4282452030868613117[5] = 0;
   out_4282452030868613117[6] = 0;
   out_4282452030868613117[7] = 0;
   out_4282452030868613117[8] = 0;
   out_4282452030868613117[9] = 0;
   out_4282452030868613117[10] = 0;
   out_4282452030868613117[11] = 0;
   out_4282452030868613117[12] = 0;
   out_4282452030868613117[13] = 1;
   out_4282452030868613117[14] = 0;
   out_4282452030868613117[15] = 0;
}
void h_30(double *state, double *unused, double *out_4340973857890209932) {
   out_4340973857890209932[0] = state[4];
}
void H_30(double *state, double *unused, double *out_9018579265503837789) {
   out_9018579265503837789[0] = 0;
   out_9018579265503837789[1] = 0;
   out_9018579265503837789[2] = 0;
   out_9018579265503837789[3] = 0;
   out_9018579265503837789[4] = 1;
   out_9018579265503837789[5] = 0;
   out_9018579265503837789[6] = 0;
   out_9018579265503837789[7] = 0;
}
void h_26(double *state, double *unused, double *out_8879862999984016737) {
   out_8879862999984016737[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5890327257682678726) {
   out_5890327257682678726[0] = 0;
   out_5890327257682678726[1] = 0;
   out_5890327257682678726[2] = 0;
   out_5890327257682678726[3] = 0;
   out_5890327257682678726[4] = 0;
   out_5890327257682678726[5] = 0;
   out_5890327257682678726[6] = 0;
   out_5890327257682678726[7] = 1;
}
void h_27(double *state, double *unused, double *out_5783911721250857558) {
   out_5783911721250857558[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7730997277667212477) {
   out_7730997277667212477[0] = 0;
   out_7730997277667212477[1] = 0;
   out_7730997277667212477[2] = 0;
   out_7730997277667212477[3] = 1;
   out_7730997277667212477[4] = 0;
   out_7730997277667212477[5] = 0;
   out_7730997277667212477[6] = 0;
   out_7730997277667212477[7] = 0;
}
void h_29(double *state, double *unused, double *out_7377759909334542617) {
   out_7377759909334542617[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1807576905715320827) {
   out_1807576905715320827[0] = 0;
   out_1807576905715320827[1] = 1;
   out_1807576905715320827[2] = 0;
   out_1807576905715320827[3] = 0;
   out_1807576905715320827[4] = 0;
   out_1807576905715320827[5] = 0;
   out_1807576905715320827[6] = 0;
   out_1807576905715320827[7] = 0;
}
void h_28(double *state, double *unused, double *out_1359095267131455886) {
   out_1359095267131455886[0] = state[5];
   out_1359095267131455886[1] = state[6];
}
void H_28(double *state, double *unused, double *out_7133737830352249877) {
   out_7133737830352249877[0] = 0;
   out_7133737830352249877[1] = 0;
   out_7133737830352249877[2] = 0;
   out_7133737830352249877[3] = 0;
   out_7133737830352249877[4] = 0;
   out_7133737830352249877[5] = 1;
   out_7133737830352249877[6] = 0;
   out_7133737830352249877[7] = 0;
   out_7133737830352249877[8] = 0;
   out_7133737830352249877[9] = 0;
   out_7133737830352249877[10] = 0;
   out_7133737830352249877[11] = 0;
   out_7133737830352249877[12] = 0;
   out_7133737830352249877[13] = 0;
   out_7133737830352249877[14] = 1;
   out_7133737830352249877[15] = 0;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_3617894249277739926) {
  err_fun(nom_x, delta_x, out_3617894249277739926);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3153453723933201675) {
  inv_err_fun(nom_x, true_x, out_3153453723933201675);
}
void car_H_mod_fun(double *state, double *out_3801013484012494542) {
  H_mod_fun(state, out_3801013484012494542);
}
void car_f_fun(double *state, double dt, double *out_7557349604050774031) {
  f_fun(state,  dt, out_7557349604050774031);
}
void car_F_fun(double *state, double dt, double *out_2200360112784492485) {
  F_fun(state,  dt, out_2200360112784492485);
}
void car_h_25(double *state, double *unused, double *out_3119584906358678393) {
  h_25(state, unused, out_3119584906358678393);
}
void car_H_25(double *state, double *unused, double *out_185734102743469453) {
  H_25(state, unused, out_185734102743469453);
}
void car_h_24(double *state, double *unused, double *out_5264692651938961980) {
  h_24(state, unused, out_5264692651938961980);
}
void car_H_24(double *state, double *unused, double *out_4282452030868613117) {
  H_24(state, unused, out_4282452030868613117);
}
void car_h_30(double *state, double *unused, double *out_4340973857890209932) {
  h_30(state, unused, out_4340973857890209932);
}
void car_H_30(double *state, double *unused, double *out_9018579265503837789) {
  H_30(state, unused, out_9018579265503837789);
}
void car_h_26(double *state, double *unused, double *out_8879862999984016737) {
  h_26(state, unused, out_8879862999984016737);
}
void car_H_26(double *state, double *unused, double *out_5890327257682678726) {
  H_26(state, unused, out_5890327257682678726);
}
void car_h_27(double *state, double *unused, double *out_5783911721250857558) {
  h_27(state, unused, out_5783911721250857558);
}
void car_H_27(double *state, double *unused, double *out_7730997277667212477) {
  H_27(state, unused, out_7730997277667212477);
}
void car_h_29(double *state, double *unused, double *out_7377759909334542617) {
  h_29(state, unused, out_7377759909334542617);
}
void car_H_29(double *state, double *unused, double *out_1807576905715320827) {
  H_29(state, unused, out_1807576905715320827);
}
void car_h_28(double *state, double *unused, double *out_1359095267131455886) {
  h_28(state, unused, out_1359095267131455886);
}
void car_H_28(double *state, double *unused, double *out_7133737830352249877) {
  H_28(state, unused, out_7133737830352249877);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
};

ekf_init(car);
