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
void err_fun(double *nom_x, double *delta_x, double *out_7817746876052121143) {
   out_7817746876052121143[0] = delta_x[0] + nom_x[0];
   out_7817746876052121143[1] = delta_x[1] + nom_x[1];
   out_7817746876052121143[2] = delta_x[2] + nom_x[2];
   out_7817746876052121143[3] = delta_x[3] + nom_x[3];
   out_7817746876052121143[4] = delta_x[4] + nom_x[4];
   out_7817746876052121143[5] = delta_x[5] + nom_x[5];
   out_7817746876052121143[6] = delta_x[6] + nom_x[6];
   out_7817746876052121143[7] = delta_x[7] + nom_x[7];
   out_7817746876052121143[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7992920226282469183) {
   out_7992920226282469183[0] = -nom_x[0] + true_x[0];
   out_7992920226282469183[1] = -nom_x[1] + true_x[1];
   out_7992920226282469183[2] = -nom_x[2] + true_x[2];
   out_7992920226282469183[3] = -nom_x[3] + true_x[3];
   out_7992920226282469183[4] = -nom_x[4] + true_x[4];
   out_7992920226282469183[5] = -nom_x[5] + true_x[5];
   out_7992920226282469183[6] = -nom_x[6] + true_x[6];
   out_7992920226282469183[7] = -nom_x[7] + true_x[7];
   out_7992920226282469183[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6026078546838174210) {
   out_6026078546838174210[0] = 1.0;
   out_6026078546838174210[1] = 0;
   out_6026078546838174210[2] = 0;
   out_6026078546838174210[3] = 0;
   out_6026078546838174210[4] = 0;
   out_6026078546838174210[5] = 0;
   out_6026078546838174210[6] = 0;
   out_6026078546838174210[7] = 0;
   out_6026078546838174210[8] = 0;
   out_6026078546838174210[9] = 0;
   out_6026078546838174210[10] = 1.0;
   out_6026078546838174210[11] = 0;
   out_6026078546838174210[12] = 0;
   out_6026078546838174210[13] = 0;
   out_6026078546838174210[14] = 0;
   out_6026078546838174210[15] = 0;
   out_6026078546838174210[16] = 0;
   out_6026078546838174210[17] = 0;
   out_6026078546838174210[18] = 0;
   out_6026078546838174210[19] = 0;
   out_6026078546838174210[20] = 1.0;
   out_6026078546838174210[21] = 0;
   out_6026078546838174210[22] = 0;
   out_6026078546838174210[23] = 0;
   out_6026078546838174210[24] = 0;
   out_6026078546838174210[25] = 0;
   out_6026078546838174210[26] = 0;
   out_6026078546838174210[27] = 0;
   out_6026078546838174210[28] = 0;
   out_6026078546838174210[29] = 0;
   out_6026078546838174210[30] = 1.0;
   out_6026078546838174210[31] = 0;
   out_6026078546838174210[32] = 0;
   out_6026078546838174210[33] = 0;
   out_6026078546838174210[34] = 0;
   out_6026078546838174210[35] = 0;
   out_6026078546838174210[36] = 0;
   out_6026078546838174210[37] = 0;
   out_6026078546838174210[38] = 0;
   out_6026078546838174210[39] = 0;
   out_6026078546838174210[40] = 1.0;
   out_6026078546838174210[41] = 0;
   out_6026078546838174210[42] = 0;
   out_6026078546838174210[43] = 0;
   out_6026078546838174210[44] = 0;
   out_6026078546838174210[45] = 0;
   out_6026078546838174210[46] = 0;
   out_6026078546838174210[47] = 0;
   out_6026078546838174210[48] = 0;
   out_6026078546838174210[49] = 0;
   out_6026078546838174210[50] = 1.0;
   out_6026078546838174210[51] = 0;
   out_6026078546838174210[52] = 0;
   out_6026078546838174210[53] = 0;
   out_6026078546838174210[54] = 0;
   out_6026078546838174210[55] = 0;
   out_6026078546838174210[56] = 0;
   out_6026078546838174210[57] = 0;
   out_6026078546838174210[58] = 0;
   out_6026078546838174210[59] = 0;
   out_6026078546838174210[60] = 1.0;
   out_6026078546838174210[61] = 0;
   out_6026078546838174210[62] = 0;
   out_6026078546838174210[63] = 0;
   out_6026078546838174210[64] = 0;
   out_6026078546838174210[65] = 0;
   out_6026078546838174210[66] = 0;
   out_6026078546838174210[67] = 0;
   out_6026078546838174210[68] = 0;
   out_6026078546838174210[69] = 0;
   out_6026078546838174210[70] = 1.0;
   out_6026078546838174210[71] = 0;
   out_6026078546838174210[72] = 0;
   out_6026078546838174210[73] = 0;
   out_6026078546838174210[74] = 0;
   out_6026078546838174210[75] = 0;
   out_6026078546838174210[76] = 0;
   out_6026078546838174210[77] = 0;
   out_6026078546838174210[78] = 0;
   out_6026078546838174210[79] = 0;
   out_6026078546838174210[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2785568340592059367) {
   out_2785568340592059367[0] = state[0];
   out_2785568340592059367[1] = state[1];
   out_2785568340592059367[2] = state[2];
   out_2785568340592059367[3] = state[3];
   out_2785568340592059367[4] = state[4];
   out_2785568340592059367[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2785568340592059367[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2785568340592059367[7] = state[7];
   out_2785568340592059367[8] = state[8];
}
void F_fun(double *state, double dt, double *out_3250476068734952827) {
   out_3250476068734952827[0] = 1;
   out_3250476068734952827[1] = 0;
   out_3250476068734952827[2] = 0;
   out_3250476068734952827[3] = 0;
   out_3250476068734952827[4] = 0;
   out_3250476068734952827[5] = 0;
   out_3250476068734952827[6] = 0;
   out_3250476068734952827[7] = 0;
   out_3250476068734952827[8] = 0;
   out_3250476068734952827[9] = 0;
   out_3250476068734952827[10] = 1;
   out_3250476068734952827[11] = 0;
   out_3250476068734952827[12] = 0;
   out_3250476068734952827[13] = 0;
   out_3250476068734952827[14] = 0;
   out_3250476068734952827[15] = 0;
   out_3250476068734952827[16] = 0;
   out_3250476068734952827[17] = 0;
   out_3250476068734952827[18] = 0;
   out_3250476068734952827[19] = 0;
   out_3250476068734952827[20] = 1;
   out_3250476068734952827[21] = 0;
   out_3250476068734952827[22] = 0;
   out_3250476068734952827[23] = 0;
   out_3250476068734952827[24] = 0;
   out_3250476068734952827[25] = 0;
   out_3250476068734952827[26] = 0;
   out_3250476068734952827[27] = 0;
   out_3250476068734952827[28] = 0;
   out_3250476068734952827[29] = 0;
   out_3250476068734952827[30] = 1;
   out_3250476068734952827[31] = 0;
   out_3250476068734952827[32] = 0;
   out_3250476068734952827[33] = 0;
   out_3250476068734952827[34] = 0;
   out_3250476068734952827[35] = 0;
   out_3250476068734952827[36] = 0;
   out_3250476068734952827[37] = 0;
   out_3250476068734952827[38] = 0;
   out_3250476068734952827[39] = 0;
   out_3250476068734952827[40] = 1;
   out_3250476068734952827[41] = 0;
   out_3250476068734952827[42] = 0;
   out_3250476068734952827[43] = 0;
   out_3250476068734952827[44] = 0;
   out_3250476068734952827[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3250476068734952827[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3250476068734952827[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3250476068734952827[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3250476068734952827[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3250476068734952827[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3250476068734952827[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3250476068734952827[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3250476068734952827[53] = -9.8000000000000007*dt;
   out_3250476068734952827[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3250476068734952827[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3250476068734952827[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3250476068734952827[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3250476068734952827[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3250476068734952827[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3250476068734952827[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3250476068734952827[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3250476068734952827[62] = 0;
   out_3250476068734952827[63] = 0;
   out_3250476068734952827[64] = 0;
   out_3250476068734952827[65] = 0;
   out_3250476068734952827[66] = 0;
   out_3250476068734952827[67] = 0;
   out_3250476068734952827[68] = 0;
   out_3250476068734952827[69] = 0;
   out_3250476068734952827[70] = 1;
   out_3250476068734952827[71] = 0;
   out_3250476068734952827[72] = 0;
   out_3250476068734952827[73] = 0;
   out_3250476068734952827[74] = 0;
   out_3250476068734952827[75] = 0;
   out_3250476068734952827[76] = 0;
   out_3250476068734952827[77] = 0;
   out_3250476068734952827[78] = 0;
   out_3250476068734952827[79] = 0;
   out_3250476068734952827[80] = 1;
}
void h_25(double *state, double *unused, double *out_7672242725304116900) {
   out_7672242725304116900[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7909775909300775539) {
   out_7909775909300775539[0] = 0;
   out_7909775909300775539[1] = 0;
   out_7909775909300775539[2] = 0;
   out_7909775909300775539[3] = 0;
   out_7909775909300775539[4] = 0;
   out_7909775909300775539[5] = 0;
   out_7909775909300775539[6] = 1;
   out_7909775909300775539[7] = 0;
   out_7909775909300775539[8] = 0;
}
void h_24(double *state, double *unused, double *out_1195134669565904685) {
   out_1195134669565904685[0] = state[4];
   out_1195134669565904685[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8364318565403276511) {
   out_8364318565403276511[0] = 0;
   out_8364318565403276511[1] = 0;
   out_8364318565403276511[2] = 0;
   out_8364318565403276511[3] = 0;
   out_8364318565403276511[4] = 1;
   out_8364318565403276511[5] = 0;
   out_8364318565403276511[6] = 0;
   out_8364318565403276511[7] = 0;
   out_8364318565403276511[8] = 0;
   out_8364318565403276511[9] = 0;
   out_8364318565403276511[10] = 0;
   out_8364318565403276511[11] = 0;
   out_8364318565403276511[12] = 0;
   out_8364318565403276511[13] = 0;
   out_8364318565403276511[14] = 1;
   out_8364318565403276511[15] = 0;
   out_8364318565403276511[16] = 0;
   out_8364318565403276511[17] = 0;
}
void h_30(double *state, double *unused, double *out_7137051820151412366) {
   out_7137051820151412366[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5391442950793526912) {
   out_5391442950793526912[0] = 0;
   out_5391442950793526912[1] = 0;
   out_5391442950793526912[2] = 0;
   out_5391442950793526912[3] = 0;
   out_5391442950793526912[4] = 1;
   out_5391442950793526912[5] = 0;
   out_5391442950793526912[6] = 0;
   out_5391442950793526912[7] = 0;
   out_5391442950793526912[8] = 0;
}
void h_26(double *state, double *unused, double *out_1844917351855428399) {
   out_1844917351855428399[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6795464845534719853) {
   out_6795464845534719853[0] = 0;
   out_6795464845534719853[1] = 0;
   out_6795464845534719853[2] = 0;
   out_6795464845534719853[3] = 0;
   out_6795464845534719853[4] = 0;
   out_6795464845534719853[5] = 0;
   out_6795464845534719853[6] = 0;
   out_6795464845534719853[7] = 1;
   out_6795464845534719853[8] = 0;
}
void h_27(double *state, double *unused, double *out_8275558609980594417) {
   out_8275558609980594417[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8232865905465111096) {
   out_8232865905465111096[0] = 0;
   out_8232865905465111096[1] = 0;
   out_8232865905465111096[2] = 0;
   out_8232865905465111096[3] = 1;
   out_8232865905465111096[4] = 0;
   out_8232865905465111096[5] = 0;
   out_8232865905465111096[6] = 0;
   out_8232865905465111096[7] = 0;
   out_8232865905465111096[8] = 0;
}
void h_29(double *state, double *unused, double *out_6533735935474934849) {
   out_6533735935474934849[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6519503178595560063) {
   out_6519503178595560063[0] = 0;
   out_6519503178595560063[1] = 1;
   out_6519503178595560063[2] = 0;
   out_6519503178595560063[3] = 0;
   out_6519503178595560063[4] = 0;
   out_6519503178595560063[5] = 0;
   out_6519503178595560063[6] = 0;
   out_6519503178595560063[7] = 0;
   out_6519503178595560063[8] = 0;
}
void h_28(double *state, double *unused, double *out_2945008839124103255) {
   out_2945008839124103255[0] = state[0];
}
void H_28(double *state, double *unused, double *out_8483133450160886314) {
   out_8483133450160886314[0] = 1;
   out_8483133450160886314[1] = 0;
   out_8483133450160886314[2] = 0;
   out_8483133450160886314[3] = 0;
   out_8483133450160886314[4] = 0;
   out_8483133450160886314[5] = 0;
   out_8483133450160886314[6] = 0;
   out_8483133450160886314[7] = 0;
   out_8483133450160886314[8] = 0;
}
void h_31(double *state, double *unused, double *out_3313942584156546391) {
   out_3313942584156546391[0] = state[8];
}
void H_31(double *state, double *unused, double *out_7879129947423815111) {
   out_7879129947423815111[0] = 0;
   out_7879129947423815111[1] = 0;
   out_7879129947423815111[2] = 0;
   out_7879129947423815111[3] = 0;
   out_7879129947423815111[4] = 0;
   out_7879129947423815111[5] = 0;
   out_7879129947423815111[6] = 0;
   out_7879129947423815111[7] = 0;
   out_7879129947423815111[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_7817746876052121143) {
  err_fun(nom_x, delta_x, out_7817746876052121143);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7992920226282469183) {
  inv_err_fun(nom_x, true_x, out_7992920226282469183);
}
void car_H_mod_fun(double *state, double *out_6026078546838174210) {
  H_mod_fun(state, out_6026078546838174210);
}
void car_f_fun(double *state, double dt, double *out_2785568340592059367) {
  f_fun(state,  dt, out_2785568340592059367);
}
void car_F_fun(double *state, double dt, double *out_3250476068734952827) {
  F_fun(state,  dt, out_3250476068734952827);
}
void car_h_25(double *state, double *unused, double *out_7672242725304116900) {
  h_25(state, unused, out_7672242725304116900);
}
void car_H_25(double *state, double *unused, double *out_7909775909300775539) {
  H_25(state, unused, out_7909775909300775539);
}
void car_h_24(double *state, double *unused, double *out_1195134669565904685) {
  h_24(state, unused, out_1195134669565904685);
}
void car_H_24(double *state, double *unused, double *out_8364318565403276511) {
  H_24(state, unused, out_8364318565403276511);
}
void car_h_30(double *state, double *unused, double *out_7137051820151412366) {
  h_30(state, unused, out_7137051820151412366);
}
void car_H_30(double *state, double *unused, double *out_5391442950793526912) {
  H_30(state, unused, out_5391442950793526912);
}
void car_h_26(double *state, double *unused, double *out_1844917351855428399) {
  h_26(state, unused, out_1844917351855428399);
}
void car_H_26(double *state, double *unused, double *out_6795464845534719853) {
  H_26(state, unused, out_6795464845534719853);
}
void car_h_27(double *state, double *unused, double *out_8275558609980594417) {
  h_27(state, unused, out_8275558609980594417);
}
void car_H_27(double *state, double *unused, double *out_8232865905465111096) {
  H_27(state, unused, out_8232865905465111096);
}
void car_h_29(double *state, double *unused, double *out_6533735935474934849) {
  h_29(state, unused, out_6533735935474934849);
}
void car_H_29(double *state, double *unused, double *out_6519503178595560063) {
  H_29(state, unused, out_6519503178595560063);
}
void car_h_28(double *state, double *unused, double *out_2945008839124103255) {
  h_28(state, unused, out_2945008839124103255);
}
void car_H_28(double *state, double *unused, double *out_8483133450160886314) {
  H_28(state, unused, out_8483133450160886314);
}
void car_h_31(double *state, double *unused, double *out_3313942584156546391) {
  h_31(state, unused, out_3313942584156546391);
}
void car_H_31(double *state, double *unused, double *out_7879129947423815111) {
  H_31(state, unused, out_7879129947423815111);
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
