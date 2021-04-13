#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_3617894249277739926);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3153453723933201675);
void car_H_mod_fun(double *state, double *out_3801013484012494542);
void car_f_fun(double *state, double dt, double *out_7557349604050774031);
void car_F_fun(double *state, double dt, double *out_2200360112784492485);
void car_h_25(double *state, double *unused, double *out_3119584906358678393);
void car_H_25(double *state, double *unused, double *out_185734102743469453);
void car_h_24(double *state, double *unused, double *out_5264692651938961980);
void car_H_24(double *state, double *unused, double *out_4282452030868613117);
void car_h_30(double *state, double *unused, double *out_4340973857890209932);
void car_H_30(double *state, double *unused, double *out_9018579265503837789);
void car_h_26(double *state, double *unused, double *out_8879862999984016737);
void car_H_26(double *state, double *unused, double *out_5890327257682678726);
void car_h_27(double *state, double *unused, double *out_5783911721250857558);
void car_H_27(double *state, double *unused, double *out_7730997277667212477);
void car_h_29(double *state, double *unused, double *out_7377759909334542617);
void car_H_29(double *state, double *unused, double *out_1807576905715320827);
void car_h_28(double *state, double *unused, double *out_1359095267131455886);
void car_H_28(double *state, double *unused, double *out_7133737830352249877);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}