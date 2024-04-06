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
void live_H(double *in_vec, double *out_4686531166020845656);
void live_err_fun(double *nom_x, double *delta_x, double *out_5460360781474494139);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_1058204392407215823);
void live_H_mod_fun(double *state, double *out_2613468595030008435);
void live_f_fun(double *state, double dt, double *out_5656405977393980005);
void live_F_fun(double *state, double dt, double *out_4808132447777396135);
void live_h_4(double *state, double *unused, double *out_6909017026966328284);
void live_H_4(double *state, double *unused, double *out_5318369041649147725);
void live_h_9(double *state, double *unused, double *out_4067103874151876147);
void live_H_9(double *state, double *unused, double *out_5841156096795956421);
void live_h_10(double *state, double *unused, double *out_878258855426047000);
void live_H_10(double *state, double *unused, double *out_4001703992801318308);
void live_h_12(double *state, double *unused, double *out_4904646255259357299);
void live_H_12(double *state, double *unused, double *out_8108918624028442096);
void live_h_35(double *state, double *unused, double *out_3108800339151075316);
void live_H_35(double *state, double *unused, double *out_8685031099021755101);
void live_h_32(double *state, double *unused, double *out_7089824626637357349);
void live_H_32(double *state, double *unused, double *out_5675977286590873415);
void live_h_13(double *state, double *unused, double *out_8589596162138588498);
void live_H_13(double *state, double *unused, double *out_1342061494263837266);
void live_h_14(double *state, double *unused, double *out_4067103874151876147);
void live_H_14(double *state, double *unused, double *out_5841156096795956421);
void live_h_33(double *state, double *unused, double *out_3949779056430118567);
void live_H_33(double *state, double *unused, double *out_6611155970048938911);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}