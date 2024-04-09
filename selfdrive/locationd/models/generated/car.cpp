#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
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
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8807790435319836399) {
   out_8807790435319836399[0] = delta_x[0] + nom_x[0];
   out_8807790435319836399[1] = delta_x[1] + nom_x[1];
   out_8807790435319836399[2] = delta_x[2] + nom_x[2];
   out_8807790435319836399[3] = delta_x[3] + nom_x[3];
   out_8807790435319836399[4] = delta_x[4] + nom_x[4];
   out_8807790435319836399[5] = delta_x[5] + nom_x[5];
   out_8807790435319836399[6] = delta_x[6] + nom_x[6];
   out_8807790435319836399[7] = delta_x[7] + nom_x[7];
   out_8807790435319836399[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2112500145296121134) {
   out_2112500145296121134[0] = -nom_x[0] + true_x[0];
   out_2112500145296121134[1] = -nom_x[1] + true_x[1];
   out_2112500145296121134[2] = -nom_x[2] + true_x[2];
   out_2112500145296121134[3] = -nom_x[3] + true_x[3];
   out_2112500145296121134[4] = -nom_x[4] + true_x[4];
   out_2112500145296121134[5] = -nom_x[5] + true_x[5];
   out_2112500145296121134[6] = -nom_x[6] + true_x[6];
   out_2112500145296121134[7] = -nom_x[7] + true_x[7];
   out_2112500145296121134[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_1713333510737573648) {
   out_1713333510737573648[0] = 1.0;
   out_1713333510737573648[1] = 0;
   out_1713333510737573648[2] = 0;
   out_1713333510737573648[3] = 0;
   out_1713333510737573648[4] = 0;
   out_1713333510737573648[5] = 0;
   out_1713333510737573648[6] = 0;
   out_1713333510737573648[7] = 0;
   out_1713333510737573648[8] = 0;
   out_1713333510737573648[9] = 0;
   out_1713333510737573648[10] = 1.0;
   out_1713333510737573648[11] = 0;
   out_1713333510737573648[12] = 0;
   out_1713333510737573648[13] = 0;
   out_1713333510737573648[14] = 0;
   out_1713333510737573648[15] = 0;
   out_1713333510737573648[16] = 0;
   out_1713333510737573648[17] = 0;
   out_1713333510737573648[18] = 0;
   out_1713333510737573648[19] = 0;
   out_1713333510737573648[20] = 1.0;
   out_1713333510737573648[21] = 0;
   out_1713333510737573648[22] = 0;
   out_1713333510737573648[23] = 0;
   out_1713333510737573648[24] = 0;
   out_1713333510737573648[25] = 0;
   out_1713333510737573648[26] = 0;
   out_1713333510737573648[27] = 0;
   out_1713333510737573648[28] = 0;
   out_1713333510737573648[29] = 0;
   out_1713333510737573648[30] = 1.0;
   out_1713333510737573648[31] = 0;
   out_1713333510737573648[32] = 0;
   out_1713333510737573648[33] = 0;
   out_1713333510737573648[34] = 0;
   out_1713333510737573648[35] = 0;
   out_1713333510737573648[36] = 0;
   out_1713333510737573648[37] = 0;
   out_1713333510737573648[38] = 0;
   out_1713333510737573648[39] = 0;
   out_1713333510737573648[40] = 1.0;
   out_1713333510737573648[41] = 0;
   out_1713333510737573648[42] = 0;
   out_1713333510737573648[43] = 0;
   out_1713333510737573648[44] = 0;
   out_1713333510737573648[45] = 0;
   out_1713333510737573648[46] = 0;
   out_1713333510737573648[47] = 0;
   out_1713333510737573648[48] = 0;
   out_1713333510737573648[49] = 0;
   out_1713333510737573648[50] = 1.0;
   out_1713333510737573648[51] = 0;
   out_1713333510737573648[52] = 0;
   out_1713333510737573648[53] = 0;
   out_1713333510737573648[54] = 0;
   out_1713333510737573648[55] = 0;
   out_1713333510737573648[56] = 0;
   out_1713333510737573648[57] = 0;
   out_1713333510737573648[58] = 0;
   out_1713333510737573648[59] = 0;
   out_1713333510737573648[60] = 1.0;
   out_1713333510737573648[61] = 0;
   out_1713333510737573648[62] = 0;
   out_1713333510737573648[63] = 0;
   out_1713333510737573648[64] = 0;
   out_1713333510737573648[65] = 0;
   out_1713333510737573648[66] = 0;
   out_1713333510737573648[67] = 0;
   out_1713333510737573648[68] = 0;
   out_1713333510737573648[69] = 0;
   out_1713333510737573648[70] = 1.0;
   out_1713333510737573648[71] = 0;
   out_1713333510737573648[72] = 0;
   out_1713333510737573648[73] = 0;
   out_1713333510737573648[74] = 0;
   out_1713333510737573648[75] = 0;
   out_1713333510737573648[76] = 0;
   out_1713333510737573648[77] = 0;
   out_1713333510737573648[78] = 0;
   out_1713333510737573648[79] = 0;
   out_1713333510737573648[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2180708203845298947) {
   out_2180708203845298947[0] = state[0];
   out_2180708203845298947[1] = state[1];
   out_2180708203845298947[2] = state[2];
   out_2180708203845298947[3] = state[3];
   out_2180708203845298947[4] = state[4];
   out_2180708203845298947[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2180708203845298947[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2180708203845298947[7] = state[7];
   out_2180708203845298947[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8567561114388139166) {
   out_8567561114388139166[0] = 1;
   out_8567561114388139166[1] = 0;
   out_8567561114388139166[2] = 0;
   out_8567561114388139166[3] = 0;
   out_8567561114388139166[4] = 0;
   out_8567561114388139166[5] = 0;
   out_8567561114388139166[6] = 0;
   out_8567561114388139166[7] = 0;
   out_8567561114388139166[8] = 0;
   out_8567561114388139166[9] = 0;
   out_8567561114388139166[10] = 1;
   out_8567561114388139166[11] = 0;
   out_8567561114388139166[12] = 0;
   out_8567561114388139166[13] = 0;
   out_8567561114388139166[14] = 0;
   out_8567561114388139166[15] = 0;
   out_8567561114388139166[16] = 0;
   out_8567561114388139166[17] = 0;
   out_8567561114388139166[18] = 0;
   out_8567561114388139166[19] = 0;
   out_8567561114388139166[20] = 1;
   out_8567561114388139166[21] = 0;
   out_8567561114388139166[22] = 0;
   out_8567561114388139166[23] = 0;
   out_8567561114388139166[24] = 0;
   out_8567561114388139166[25] = 0;
   out_8567561114388139166[26] = 0;
   out_8567561114388139166[27] = 0;
   out_8567561114388139166[28] = 0;
   out_8567561114388139166[29] = 0;
   out_8567561114388139166[30] = 1;
   out_8567561114388139166[31] = 0;
   out_8567561114388139166[32] = 0;
   out_8567561114388139166[33] = 0;
   out_8567561114388139166[34] = 0;
   out_8567561114388139166[35] = 0;
   out_8567561114388139166[36] = 0;
   out_8567561114388139166[37] = 0;
   out_8567561114388139166[38] = 0;
   out_8567561114388139166[39] = 0;
   out_8567561114388139166[40] = 1;
   out_8567561114388139166[41] = 0;
   out_8567561114388139166[42] = 0;
   out_8567561114388139166[43] = 0;
   out_8567561114388139166[44] = 0;
   out_8567561114388139166[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8567561114388139166[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8567561114388139166[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8567561114388139166[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8567561114388139166[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8567561114388139166[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8567561114388139166[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8567561114388139166[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8567561114388139166[53] = -9.8000000000000007*dt;
   out_8567561114388139166[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8567561114388139166[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8567561114388139166[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8567561114388139166[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8567561114388139166[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8567561114388139166[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8567561114388139166[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8567561114388139166[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8567561114388139166[62] = 0;
   out_8567561114388139166[63] = 0;
   out_8567561114388139166[64] = 0;
   out_8567561114388139166[65] = 0;
   out_8567561114388139166[66] = 0;
   out_8567561114388139166[67] = 0;
   out_8567561114388139166[68] = 0;
   out_8567561114388139166[69] = 0;
   out_8567561114388139166[70] = 1;
   out_8567561114388139166[71] = 0;
   out_8567561114388139166[72] = 0;
   out_8567561114388139166[73] = 0;
   out_8567561114388139166[74] = 0;
   out_8567561114388139166[75] = 0;
   out_8567561114388139166[76] = 0;
   out_8567561114388139166[77] = 0;
   out_8567561114388139166[78] = 0;
   out_8567561114388139166[79] = 0;
   out_8567561114388139166[80] = 1;
}
void h_25(double *state, double *unused, double *out_7592369920789948917) {
   out_7592369920789948917[0] = state[6];
}
void H_25(double *state, double *unused, double *out_9130032943935932027) {
   out_9130032943935932027[0] = 0;
   out_9130032943935932027[1] = 0;
   out_9130032943935932027[2] = 0;
   out_9130032943935932027[3] = 0;
   out_9130032943935932027[4] = 0;
   out_9130032943935932027[5] = 0;
   out_9130032943935932027[6] = 1;
   out_9130032943935932027[7] = 0;
   out_9130032943935932027[8] = 0;
}
void h_24(double *state, double *unused, double *out_7234279162208326407) {
   out_7234279162208326407[0] = state[4];
   out_7234279162208326407[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6957383344930432461) {
   out_6957383344930432461[0] = 0;
   out_6957383344930432461[1] = 0;
   out_6957383344930432461[2] = 0;
   out_6957383344930432461[3] = 0;
   out_6957383344930432461[4] = 1;
   out_6957383344930432461[5] = 0;
   out_6957383344930432461[6] = 0;
   out_6957383344930432461[7] = 0;
   out_6957383344930432461[8] = 0;
   out_6957383344930432461[9] = 0;
   out_6957383344930432461[10] = 0;
   out_6957383344930432461[11] = 0;
   out_6957383344930432461[12] = 0;
   out_6957383344930432461[13] = 0;
   out_6957383344930432461[14] = 1;
   out_6957383344930432461[15] = 0;
   out_6957383344930432461[16] = 0;
   out_6957383344930432461[17] = 0;
}
void h_30(double *state, double *unused, double *out_7867563983074454806) {
   out_7867563983074454806[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6798378171266370962) {
   out_6798378171266370962[0] = 0;
   out_6798378171266370962[1] = 0;
   out_6798378171266370962[2] = 0;
   out_6798378171266370962[3] = 0;
   out_6798378171266370962[4] = 1;
   out_6798378171266370962[5] = 0;
   out_6798378171266370962[6] = 0;
   out_6798378171266370962[7] = 0;
   out_6798378171266370962[8] = 0;
}
void h_26(double *state, double *unused, double *out_2865890715413609161) {
   out_2865890715413609161[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5388529625061875803) {
   out_5388529625061875803[0] = 0;
   out_5388529625061875803[1] = 0;
   out_5388529625061875803[2] = 0;
   out_5388529625061875803[3] = 0;
   out_5388529625061875803[4] = 0;
   out_5388529625061875803[5] = 0;
   out_5388529625061875803[6] = 0;
   out_5388529625061875803[7] = 1;
   out_5388529625061875803[8] = 0;
}
void h_27(double *state, double *unused, double *out_589607216996883138) {
   out_589607216996883138[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4574784100082427745) {
   out_4574784100082427745[0] = 0;
   out_4574784100082427745[1] = 0;
   out_4574784100082427745[2] = 0;
   out_4574784100082427745[3] = 1;
   out_4574784100082427745[4] = 0;
   out_4574784100082427745[5] = 0;
   out_4574784100082427745[6] = 0;
   out_4574784100082427745[7] = 0;
   out_4574784100082427745[8] = 0;
}
void h_29(double *state, double *unused, double *out_4936450219554139975) {
   out_4936450219554139975[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6288146826951978778) {
   out_6288146826951978778[0] = 0;
   out_6288146826951978778[1] = 1;
   out_6288146826951978778[2] = 0;
   out_6288146826951978778[3] = 0;
   out_6288146826951978778[4] = 0;
   out_6288146826951978778[5] = 0;
   out_6288146826951978778[6] = 0;
   out_6288146826951978778[7] = 0;
   out_6288146826951978778[8] = 0;
}
void h_28(double *state, double *unused, double *out_353902635904779297) {
   out_353902635904779297[0] = state[0];
}
void H_28(double *state, double *unused, double *out_7076198229688042264) {
   out_7076198229688042264[0] = 1;
   out_7076198229688042264[1] = 0;
   out_7076198229688042264[2] = 0;
   out_7076198229688042264[3] = 0;
   out_7076198229688042264[4] = 0;
   out_7076198229688042264[5] = 0;
   out_7076198229688042264[6] = 0;
   out_7076198229688042264[7] = 0;
   out_7076198229688042264[8] = 0;
}
void h_31(double *state, double *unused, double *out_2961434658035079340) {
   out_2961434658035079340[0] = state[8];
}
void H_31(double *state, double *unused, double *out_9160678905812892455) {
   out_9160678905812892455[0] = 0;
   out_9160678905812892455[1] = 0;
   out_9160678905812892455[2] = 0;
   out_9160678905812892455[3] = 0;
   out_9160678905812892455[4] = 0;
   out_9160678905812892455[5] = 0;
   out_9160678905812892455[6] = 0;
   out_9160678905812892455[7] = 0;
   out_9160678905812892455[8] = 1;
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
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_8807790435319836399) {
  err_fun(nom_x, delta_x, out_8807790435319836399);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2112500145296121134) {
  inv_err_fun(nom_x, true_x, out_2112500145296121134);
}
void car_H_mod_fun(double *state, double *out_1713333510737573648) {
  H_mod_fun(state, out_1713333510737573648);
}
void car_f_fun(double *state, double dt, double *out_2180708203845298947) {
  f_fun(state,  dt, out_2180708203845298947);
}
void car_F_fun(double *state, double dt, double *out_8567561114388139166) {
  F_fun(state,  dt, out_8567561114388139166);
}
void car_h_25(double *state, double *unused, double *out_7592369920789948917) {
  h_25(state, unused, out_7592369920789948917);
}
void car_H_25(double *state, double *unused, double *out_9130032943935932027) {
  H_25(state, unused, out_9130032943935932027);
}
void car_h_24(double *state, double *unused, double *out_7234279162208326407) {
  h_24(state, unused, out_7234279162208326407);
}
void car_H_24(double *state, double *unused, double *out_6957383344930432461) {
  H_24(state, unused, out_6957383344930432461);
}
void car_h_30(double *state, double *unused, double *out_7867563983074454806) {
  h_30(state, unused, out_7867563983074454806);
}
void car_H_30(double *state, double *unused, double *out_6798378171266370962) {
  H_30(state, unused, out_6798378171266370962);
}
void car_h_26(double *state, double *unused, double *out_2865890715413609161) {
  h_26(state, unused, out_2865890715413609161);
}
void car_H_26(double *state, double *unused, double *out_5388529625061875803) {
  H_26(state, unused, out_5388529625061875803);
}
void car_h_27(double *state, double *unused, double *out_589607216996883138) {
  h_27(state, unused, out_589607216996883138);
}
void car_H_27(double *state, double *unused, double *out_4574784100082427745) {
  H_27(state, unused, out_4574784100082427745);
}
void car_h_29(double *state, double *unused, double *out_4936450219554139975) {
  h_29(state, unused, out_4936450219554139975);
}
void car_H_29(double *state, double *unused, double *out_6288146826951978778) {
  H_29(state, unused, out_6288146826951978778);
}
void car_h_28(double *state, double *unused, double *out_353902635904779297) {
  h_28(state, unused, out_353902635904779297);
}
void car_H_28(double *state, double *unused, double *out_7076198229688042264) {
  H_28(state, unused, out_7076198229688042264);
}
void car_h_31(double *state, double *unused, double *out_2961434658035079340) {
  h_31(state, unused, out_2961434658035079340);
}
void car_H_31(double *state, double *unused, double *out_9160678905812892455) {
  H_31(state, unused, out_9160678905812892455);
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
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
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
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
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
  .extra_routines = {
  },
};

ekf_lib_init(car)
