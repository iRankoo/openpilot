#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_8807790435319836399);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2112500145296121134);
void car_H_mod_fun(double *state, double *out_1713333510737573648);
void car_f_fun(double *state, double dt, double *out_2180708203845298947);
void car_F_fun(double *state, double dt, double *out_8567561114388139166);
void car_h_25(double *state, double *unused, double *out_7592369920789948917);
void car_H_25(double *state, double *unused, double *out_9130032943935932027);
void car_h_24(double *state, double *unused, double *out_7234279162208326407);
void car_H_24(double *state, double *unused, double *out_6957383344930432461);
void car_h_30(double *state, double *unused, double *out_7867563983074454806);
void car_H_30(double *state, double *unused, double *out_6798378171266370962);
void car_h_26(double *state, double *unused, double *out_2865890715413609161);
void car_H_26(double *state, double *unused, double *out_5388529625061875803);
void car_h_27(double *state, double *unused, double *out_589607216996883138);
void car_H_27(double *state, double *unused, double *out_4574784100082427745);
void car_h_29(double *state, double *unused, double *out_4936450219554139975);
void car_H_29(double *state, double *unused, double *out_6288146826951978778);
void car_h_28(double *state, double *unused, double *out_353902635904779297);
void car_H_28(double *state, double *unused, double *out_7076198229688042264);
void car_h_31(double *state, double *unused, double *out_2961434658035079340);
void car_H_31(double *state, double *unused, double *out_9160678905812892455);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}