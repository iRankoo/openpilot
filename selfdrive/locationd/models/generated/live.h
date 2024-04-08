#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_7067658687324008333);
void live_err_fun(double *nom_x, double *delta_x, double *out_4190484818694994195);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_4666604619517225884);
void live_H_mod_fun(double *state, double *out_832045067153363880);
void live_f_fun(double *state, double dt, double *out_2532795489537895138);
void live_F_fun(double *state, double dt, double *out_6201617391151166938);
void live_h_4(double *state, double *unused, double *out_4076930992399196000);
void live_H_4(double *state, double *unused, double *out_2130574254833149188);
void live_h_9(double *state, double *unused, double *out_6916800482633854831);
void live_H_9(double *state, double *unused, double *out_2371763901462739833);
void live_h_10(double *state, double *unused, double *out_7220059633906934465);
void live_H_10(double *state, double *unused, double *out_7717948348674522949);
void live_h_12(double *state, double *unused, double *out_3944459411321133777);
void live_H_12(double *state, double *unused, double *out_7150030662865110983);
void live_h_35(double *state, double *unused, double *out_6809259272958510517);
void live_H_35(double *state, double *unused, double *out_5497236312205756564);
void live_h_32(double *state, double *unused, double *out_7709909724760785904);
void live_H_32(double *state, double *unused, double *out_1417628191488359792);
void live_h_13(double *state, double *unused, double *out_2899346340352401221);
void live_H_13(double *state, double *unused, double *out_4869931875171865865);
void live_h_14(double *state, double *unused, double *out_6916800482633854831);
void live_H_14(double *state, double *unused, double *out_2371763901462739833);
void live_h_33(double *state, double *unused, double *out_8064077977877237954);
void live_H_33(double *state, double *unused, double *out_8647793316844614168);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}