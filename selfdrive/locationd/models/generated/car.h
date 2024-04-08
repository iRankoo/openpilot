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
void car_err_fun(double *nom_x, double *delta_x, double *out_2055074247431952110);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_710933572915437639);
void car_H_mod_fun(double *state, double *out_1065691078985214700);
void car_f_fun(double *state, double dt, double *out_51088722969000182);
void car_F_fun(double *state, double dt, double *out_8609222397993749237);
void car_h_25(double *state, double *unused, double *out_6447545910113066622);
void car_H_25(double *state, double *unused, double *out_565198808842467904);
void car_h_24(double *state, double *unused, double *out_2189217925543120910);
void car_H_24(double *state, double *unused, double *out_418867061998845764);
void car_h_30(double *state, double *unused, double *out_5045581451962609654);
void car_H_30(double *state, double *unused, double *out_6351491532649148851);
void car_h_26(double *state, double *unused, double *out_7809466620754036838);
void car_H_26(double *state, double *unused, double *out_4306702127716524128);
void car_h_27(double *state, double *unused, double *out_240388122632445567);
void car_H_27(double *state, double *unused, double *out_4176728220848723940);
void car_h_29(double *state, double *unused, double *out_8256272558631294425);
void car_H_29(double *state, double *unused, double *out_6861722876963541035);
void car_h_28(double *state, double *unused, double *out_1222472215967743679);
void car_H_28(double *state, double *unused, double *out_2619033523090357667);
void car_h_31(double *state, double *unused, double *out_5036479207312905967);
void car_H_31(double *state, double *unused, double *out_534552846965507476);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}