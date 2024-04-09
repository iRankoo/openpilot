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
void live_H(double *in_vec, double *out_5284371152158479112);
void live_err_fun(double *nom_x, double *delta_x, double *out_8878190256849491120);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_6398751028731537898);
void live_H_mod_fun(double *state, double *out_2043043155980843567);
void live_f_fun(double *state, double dt, double *out_2010583831986639410);
void live_F_fun(double *state, double dt, double *out_6122016632954703737);
void live_h_4(double *state, double *unused, double *out_9003365468988045872);
void live_H_4(double *state, double *unused, double *out_875493743705787679);
void live_h_9(double *state, double *unused, double *out_7294275750039166322);
void live_H_9(double *state, double *unused, double *out_6411725191558659791);
void live_h_10(double *state, double *unused, double *out_8289912668181558020);
void live_H_10(double *state, double *unused, double *out_4756992804745523366);
void live_h_12(double *state, double *unused, double *out_2787968462159155653);
void live_H_12(double *state, double *unused, double *out_4143962664326174116);
void live_h_35(double *state, double *unused, double *out_1938642465562957644);
void live_H_35(double *state, double *unused, double *out_2491168313666819697);
void live_h_32(double *state, double *unused, double *out_6691864520126804195);
void live_H_32(double *state, double *unused, double *out_395313109181634799);
void live_h_13(double *state, double *unused, double *out_3837576363831715836);
void live_H_13(double *state, double *unused, double *out_3742987210220696894);
void live_h_14(double *state, double *unused, double *out_7294275750039166322);
void live_H_14(double *state, double *unused, double *out_6411725191558659791);
void live_h_33(double *state, double *unused, double *out_7556352649106866925);
void live_H_33(double *state, double *unused, double *out_5641725318305677301);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}