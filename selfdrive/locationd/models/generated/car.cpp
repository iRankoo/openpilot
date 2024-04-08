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
void err_fun(double *nom_x, double *delta_x, double *out_2055074247431952110) {
   out_2055074247431952110[0] = delta_x[0] + nom_x[0];
   out_2055074247431952110[1] = delta_x[1] + nom_x[1];
   out_2055074247431952110[2] = delta_x[2] + nom_x[2];
   out_2055074247431952110[3] = delta_x[3] + nom_x[3];
   out_2055074247431952110[4] = delta_x[4] + nom_x[4];
   out_2055074247431952110[5] = delta_x[5] + nom_x[5];
   out_2055074247431952110[6] = delta_x[6] + nom_x[6];
   out_2055074247431952110[7] = delta_x[7] + nom_x[7];
   out_2055074247431952110[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_710933572915437639) {
   out_710933572915437639[0] = -nom_x[0] + true_x[0];
   out_710933572915437639[1] = -nom_x[1] + true_x[1];
   out_710933572915437639[2] = -nom_x[2] + true_x[2];
   out_710933572915437639[3] = -nom_x[3] + true_x[3];
   out_710933572915437639[4] = -nom_x[4] + true_x[4];
   out_710933572915437639[5] = -nom_x[5] + true_x[5];
   out_710933572915437639[6] = -nom_x[6] + true_x[6];
   out_710933572915437639[7] = -nom_x[7] + true_x[7];
   out_710933572915437639[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_1065691078985214700) {
   out_1065691078985214700[0] = 1.0;
   out_1065691078985214700[1] = 0;
   out_1065691078985214700[2] = 0;
   out_1065691078985214700[3] = 0;
   out_1065691078985214700[4] = 0;
   out_1065691078985214700[5] = 0;
   out_1065691078985214700[6] = 0;
   out_1065691078985214700[7] = 0;
   out_1065691078985214700[8] = 0;
   out_1065691078985214700[9] = 0;
   out_1065691078985214700[10] = 1.0;
   out_1065691078985214700[11] = 0;
   out_1065691078985214700[12] = 0;
   out_1065691078985214700[13] = 0;
   out_1065691078985214700[14] = 0;
   out_1065691078985214700[15] = 0;
   out_1065691078985214700[16] = 0;
   out_1065691078985214700[17] = 0;
   out_1065691078985214700[18] = 0;
   out_1065691078985214700[19] = 0;
   out_1065691078985214700[20] = 1.0;
   out_1065691078985214700[21] = 0;
   out_1065691078985214700[22] = 0;
   out_1065691078985214700[23] = 0;
   out_1065691078985214700[24] = 0;
   out_1065691078985214700[25] = 0;
   out_1065691078985214700[26] = 0;
   out_1065691078985214700[27] = 0;
   out_1065691078985214700[28] = 0;
   out_1065691078985214700[29] = 0;
   out_1065691078985214700[30] = 1.0;
   out_1065691078985214700[31] = 0;
   out_1065691078985214700[32] = 0;
   out_1065691078985214700[33] = 0;
   out_1065691078985214700[34] = 0;
   out_1065691078985214700[35] = 0;
   out_1065691078985214700[36] = 0;
   out_1065691078985214700[37] = 0;
   out_1065691078985214700[38] = 0;
   out_1065691078985214700[39] = 0;
   out_1065691078985214700[40] = 1.0;
   out_1065691078985214700[41] = 0;
   out_1065691078985214700[42] = 0;
   out_1065691078985214700[43] = 0;
   out_1065691078985214700[44] = 0;
   out_1065691078985214700[45] = 0;
   out_1065691078985214700[46] = 0;
   out_1065691078985214700[47] = 0;
   out_1065691078985214700[48] = 0;
   out_1065691078985214700[49] = 0;
   out_1065691078985214700[50] = 1.0;
   out_1065691078985214700[51] = 0;
   out_1065691078985214700[52] = 0;
   out_1065691078985214700[53] = 0;
   out_1065691078985214700[54] = 0;
   out_1065691078985214700[55] = 0;
   out_1065691078985214700[56] = 0;
   out_1065691078985214700[57] = 0;
   out_1065691078985214700[58] = 0;
   out_1065691078985214700[59] = 0;
   out_1065691078985214700[60] = 1.0;
   out_1065691078985214700[61] = 0;
   out_1065691078985214700[62] = 0;
   out_1065691078985214700[63] = 0;
   out_1065691078985214700[64] = 0;
   out_1065691078985214700[65] = 0;
   out_1065691078985214700[66] = 0;
   out_1065691078985214700[67] = 0;
   out_1065691078985214700[68] = 0;
   out_1065691078985214700[69] = 0;
   out_1065691078985214700[70] = 1.0;
   out_1065691078985214700[71] = 0;
   out_1065691078985214700[72] = 0;
   out_1065691078985214700[73] = 0;
   out_1065691078985214700[74] = 0;
   out_1065691078985214700[75] = 0;
   out_1065691078985214700[76] = 0;
   out_1065691078985214700[77] = 0;
   out_1065691078985214700[78] = 0;
   out_1065691078985214700[79] = 0;
   out_1065691078985214700[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_51088722969000182) {
   out_51088722969000182[0] = state[0];
   out_51088722969000182[1] = state[1];
   out_51088722969000182[2] = state[2];
   out_51088722969000182[3] = state[3];
   out_51088722969000182[4] = state[4];
   out_51088722969000182[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_51088722969000182[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_51088722969000182[7] = state[7];
   out_51088722969000182[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8609222397993749237) {
   out_8609222397993749237[0] = 1;
   out_8609222397993749237[1] = 0;
   out_8609222397993749237[2] = 0;
   out_8609222397993749237[3] = 0;
   out_8609222397993749237[4] = 0;
   out_8609222397993749237[5] = 0;
   out_8609222397993749237[6] = 0;
   out_8609222397993749237[7] = 0;
   out_8609222397993749237[8] = 0;
   out_8609222397993749237[9] = 0;
   out_8609222397993749237[10] = 1;
   out_8609222397993749237[11] = 0;
   out_8609222397993749237[12] = 0;
   out_8609222397993749237[13] = 0;
   out_8609222397993749237[14] = 0;
   out_8609222397993749237[15] = 0;
   out_8609222397993749237[16] = 0;
   out_8609222397993749237[17] = 0;
   out_8609222397993749237[18] = 0;
   out_8609222397993749237[19] = 0;
   out_8609222397993749237[20] = 1;
   out_8609222397993749237[21] = 0;
   out_8609222397993749237[22] = 0;
   out_8609222397993749237[23] = 0;
   out_8609222397993749237[24] = 0;
   out_8609222397993749237[25] = 0;
   out_8609222397993749237[26] = 0;
   out_8609222397993749237[27] = 0;
   out_8609222397993749237[28] = 0;
   out_8609222397993749237[29] = 0;
   out_8609222397993749237[30] = 1;
   out_8609222397993749237[31] = 0;
   out_8609222397993749237[32] = 0;
   out_8609222397993749237[33] = 0;
   out_8609222397993749237[34] = 0;
   out_8609222397993749237[35] = 0;
   out_8609222397993749237[36] = 0;
   out_8609222397993749237[37] = 0;
   out_8609222397993749237[38] = 0;
   out_8609222397993749237[39] = 0;
   out_8609222397993749237[40] = 1;
   out_8609222397993749237[41] = 0;
   out_8609222397993749237[42] = 0;
   out_8609222397993749237[43] = 0;
   out_8609222397993749237[44] = 0;
   out_8609222397993749237[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8609222397993749237[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8609222397993749237[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8609222397993749237[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8609222397993749237[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8609222397993749237[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8609222397993749237[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8609222397993749237[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8609222397993749237[53] = -9.8000000000000007*dt;
   out_8609222397993749237[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8609222397993749237[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8609222397993749237[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8609222397993749237[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8609222397993749237[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8609222397993749237[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8609222397993749237[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8609222397993749237[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8609222397993749237[62] = 0;
   out_8609222397993749237[63] = 0;
   out_8609222397993749237[64] = 0;
   out_8609222397993749237[65] = 0;
   out_8609222397993749237[66] = 0;
   out_8609222397993749237[67] = 0;
   out_8609222397993749237[68] = 0;
   out_8609222397993749237[69] = 0;
   out_8609222397993749237[70] = 1;
   out_8609222397993749237[71] = 0;
   out_8609222397993749237[72] = 0;
   out_8609222397993749237[73] = 0;
   out_8609222397993749237[74] = 0;
   out_8609222397993749237[75] = 0;
   out_8609222397993749237[76] = 0;
   out_8609222397993749237[77] = 0;
   out_8609222397993749237[78] = 0;
   out_8609222397993749237[79] = 0;
   out_8609222397993749237[80] = 1;
}
void h_25(double *state, double *unused, double *out_6447545910113066622) {
   out_6447545910113066622[0] = state[6];
}
void H_25(double *state, double *unused, double *out_565198808842467904) {
   out_565198808842467904[0] = 0;
   out_565198808842467904[1] = 0;
   out_565198808842467904[2] = 0;
   out_565198808842467904[3] = 0;
   out_565198808842467904[4] = 0;
   out_565198808842467904[5] = 0;
   out_565198808842467904[6] = 1;
   out_565198808842467904[7] = 0;
   out_565198808842467904[8] = 0;
}
void h_24(double *state, double *unused, double *out_2189217925543120910) {
   out_2189217925543120910[0] = state[4];
   out_2189217925543120910[1] = state[5];
}
void H_24(double *state, double *unused, double *out_418867061998845764) {
   out_418867061998845764[0] = 0;
   out_418867061998845764[1] = 0;
   out_418867061998845764[2] = 0;
   out_418867061998845764[3] = 0;
   out_418867061998845764[4] = 1;
   out_418867061998845764[5] = 0;
   out_418867061998845764[6] = 0;
   out_418867061998845764[7] = 0;
   out_418867061998845764[8] = 0;
   out_418867061998845764[9] = 0;
   out_418867061998845764[10] = 0;
   out_418867061998845764[11] = 0;
   out_418867061998845764[12] = 0;
   out_418867061998845764[13] = 0;
   out_418867061998845764[14] = 1;
   out_418867061998845764[15] = 0;
   out_418867061998845764[16] = 0;
   out_418867061998845764[17] = 0;
}
void h_30(double *state, double *unused, double *out_5045581451962609654) {
   out_5045581451962609654[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6351491532649148851) {
   out_6351491532649148851[0] = 0;
   out_6351491532649148851[1] = 0;
   out_6351491532649148851[2] = 0;
   out_6351491532649148851[3] = 0;
   out_6351491532649148851[4] = 1;
   out_6351491532649148851[5] = 0;
   out_6351491532649148851[6] = 0;
   out_6351491532649148851[7] = 0;
   out_6351491532649148851[8] = 0;
}
void h_26(double *state, double *unused, double *out_7809466620754036838) {
   out_7809466620754036838[0] = state[7];
}
void H_26(double *state, double *unused, double *out_4306702127716524128) {
   out_4306702127716524128[0] = 0;
   out_4306702127716524128[1] = 0;
   out_4306702127716524128[2] = 0;
   out_4306702127716524128[3] = 0;
   out_4306702127716524128[4] = 0;
   out_4306702127716524128[5] = 0;
   out_4306702127716524128[6] = 0;
   out_4306702127716524128[7] = 1;
   out_4306702127716524128[8] = 0;
}
void h_27(double *state, double *unused, double *out_240388122632445567) {
   out_240388122632445567[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4176728220848723940) {
   out_4176728220848723940[0] = 0;
   out_4176728220848723940[1] = 0;
   out_4176728220848723940[2] = 0;
   out_4176728220848723940[3] = 1;
   out_4176728220848723940[4] = 0;
   out_4176728220848723940[5] = 0;
   out_4176728220848723940[6] = 0;
   out_4176728220848723940[7] = 0;
   out_4176728220848723940[8] = 0;
}
void h_29(double *state, double *unused, double *out_8256272558631294425) {
   out_8256272558631294425[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6861722876963541035) {
   out_6861722876963541035[0] = 0;
   out_6861722876963541035[1] = 1;
   out_6861722876963541035[2] = 0;
   out_6861722876963541035[3] = 0;
   out_6861722876963541035[4] = 0;
   out_6861722876963541035[5] = 0;
   out_6861722876963541035[6] = 0;
   out_6861722876963541035[7] = 0;
   out_6861722876963541035[8] = 0;
}
void h_28(double *state, double *unused, double *out_1222472215967743679) {
   out_1222472215967743679[0] = state[0];
}
void H_28(double *state, double *unused, double *out_2619033523090357667) {
   out_2619033523090357667[0] = 1;
   out_2619033523090357667[1] = 0;
   out_2619033523090357667[2] = 0;
   out_2619033523090357667[3] = 0;
   out_2619033523090357667[4] = 0;
   out_2619033523090357667[5] = 0;
   out_2619033523090357667[6] = 0;
   out_2619033523090357667[7] = 0;
   out_2619033523090357667[8] = 0;
}
void h_31(double *state, double *unused, double *out_5036479207312905967) {
   out_5036479207312905967[0] = state[8];
}
void H_31(double *state, double *unused, double *out_534552846965507476) {
   out_534552846965507476[0] = 0;
   out_534552846965507476[1] = 0;
   out_534552846965507476[2] = 0;
   out_534552846965507476[3] = 0;
   out_534552846965507476[4] = 0;
   out_534552846965507476[5] = 0;
   out_534552846965507476[6] = 0;
   out_534552846965507476[7] = 0;
   out_534552846965507476[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_2055074247431952110) {
  err_fun(nom_x, delta_x, out_2055074247431952110);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_710933572915437639) {
  inv_err_fun(nom_x, true_x, out_710933572915437639);
}
void car_H_mod_fun(double *state, double *out_1065691078985214700) {
  H_mod_fun(state, out_1065691078985214700);
}
void car_f_fun(double *state, double dt, double *out_51088722969000182) {
  f_fun(state,  dt, out_51088722969000182);
}
void car_F_fun(double *state, double dt, double *out_8609222397993749237) {
  F_fun(state,  dt, out_8609222397993749237);
}
void car_h_25(double *state, double *unused, double *out_6447545910113066622) {
  h_25(state, unused, out_6447545910113066622);
}
void car_H_25(double *state, double *unused, double *out_565198808842467904) {
  H_25(state, unused, out_565198808842467904);
}
void car_h_24(double *state, double *unused, double *out_2189217925543120910) {
  h_24(state, unused, out_2189217925543120910);
}
void car_H_24(double *state, double *unused, double *out_418867061998845764) {
  H_24(state, unused, out_418867061998845764);
}
void car_h_30(double *state, double *unused, double *out_5045581451962609654) {
  h_30(state, unused, out_5045581451962609654);
}
void car_H_30(double *state, double *unused, double *out_6351491532649148851) {
  H_30(state, unused, out_6351491532649148851);
}
void car_h_26(double *state, double *unused, double *out_7809466620754036838) {
  h_26(state, unused, out_7809466620754036838);
}
void car_H_26(double *state, double *unused, double *out_4306702127716524128) {
  H_26(state, unused, out_4306702127716524128);
}
void car_h_27(double *state, double *unused, double *out_240388122632445567) {
  h_27(state, unused, out_240388122632445567);
}
void car_H_27(double *state, double *unused, double *out_4176728220848723940) {
  H_27(state, unused, out_4176728220848723940);
}
void car_h_29(double *state, double *unused, double *out_8256272558631294425) {
  h_29(state, unused, out_8256272558631294425);
}
void car_H_29(double *state, double *unused, double *out_6861722876963541035) {
  H_29(state, unused, out_6861722876963541035);
}
void car_h_28(double *state, double *unused, double *out_1222472215967743679) {
  h_28(state, unused, out_1222472215967743679);
}
void car_H_28(double *state, double *unused, double *out_2619033523090357667) {
  H_28(state, unused, out_2619033523090357667);
}
void car_h_31(double *state, double *unused, double *out_5036479207312905967) {
  h_31(state, unused, out_5036479207312905967);
}
void car_H_31(double *state, double *unused, double *out_534552846965507476) {
  H_31(state, unused, out_534552846965507476);
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
