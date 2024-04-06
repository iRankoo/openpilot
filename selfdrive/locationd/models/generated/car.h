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
void car_err_fun(double *nom_x, double *delta_x, double *out_7817746876052121143);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7992920226282469183);
void car_H_mod_fun(double *state, double *out_6026078546838174210);
void car_f_fun(double *state, double dt, double *out_2785568340592059367);
void car_F_fun(double *state, double dt, double *out_3250476068734952827);
void car_h_25(double *state, double *unused, double *out_7672242725304116900);
void car_H_25(double *state, double *unused, double *out_7909775909300775539);
void car_h_24(double *state, double *unused, double *out_1195134669565904685);
void car_H_24(double *state, double *unused, double *out_8364318565403276511);
void car_h_30(double *state, double *unused, double *out_7137051820151412366);
void car_H_30(double *state, double *unused, double *out_5391442950793526912);
void car_h_26(double *state, double *unused, double *out_1844917351855428399);
void car_H_26(double *state, double *unused, double *out_6795464845534719853);
void car_h_27(double *state, double *unused, double *out_8275558609980594417);
void car_H_27(double *state, double *unused, double *out_8232865905465111096);
void car_h_29(double *state, double *unused, double *out_6533735935474934849);
void car_H_29(double *state, double *unused, double *out_6519503178595560063);
void car_h_28(double *state, double *unused, double *out_2945008839124103255);
void car_H_28(double *state, double *unused, double *out_8483133450160886314);
void car_h_31(double *state, double *unused, double *out_3313942584156546391);
void car_H_31(double *state, double *unused, double *out_7879129947423815111);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}