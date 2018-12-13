/*
 * dual_quadratic_interpolation.h
 *
 *  Created on: Dec 6, 2018
 *      Author: winona.wei
 */
#ifndef DUAL_QUADRATIC_INTERPOLATION_H
#define DUAL_QUADRATIC_INTERPOLATION_H

void dual_quadratic_interpolation(double in_x[],int xn,double in_y[],int yn,double z[6][6],double h,double t,double x,double y,double *out_p);
#endif
