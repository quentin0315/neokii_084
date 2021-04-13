#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_3(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_19(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_err_fun(double *nom_x, double *delta_x, double *out_2996855827947511820);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_6690315078635622056);
void live_H_mod_fun(double *state, double *out_3694952956564961820);
void live_f_fun(double *state, double dt, double *out_454423987521151353);
void live_F_fun(double *state, double dt, double *out_6126415230543029209);
void live_h_3(double *state, double *unused, double *out_8746128282463051673);
void live_H_3(double *state, double *unused, double *out_1720070419806077995);
void live_h_4(double *state, double *unused, double *out_7393778967491304424);
void live_H_4(double *state, double *unused, double *out_5292668717299202330);
void live_h_9(double *state, double *unused, double *out_5474774614289376666);
void live_H_9(double *state, double *unused, double *out_7615174261170058300);
void live_h_10(double *state, double *unused, double *out_7698823489387234420);
void live_H_10(double *state, double *unused, double *out_113991217274959440);
void live_h_12(double *state, double *unused, double *out_4582374479334988964);
void live_H_12(double *state, double *unused, double *out_7399224229038411899);
void live_h_31(double *state, double *unused, double *out_5103758523791028526);
void live_H_31(double *state, double *unused, double *out_4724134453894992078);
void live_h_32(double *state, double *unused, double *out_7179083296985259106);
void live_H_32(double *state, double *unused, double *out_5869912591553103052);
void live_h_13(double *state, double *unused, double *out_3457591238261864803);
void live_H_13(double *state, double *unused, double *out_9086417574309775559);
void live_h_14(double *state, double *unused, double *out_5474774614289376666);
void live_H_14(double *state, double *unused, double *out_7615174261170058300);
void live_h_19(double *state, double *unused, double *out_918413708724762745);
void live_H_19(double *state, double *unused, double *out_6398889961437062135);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}