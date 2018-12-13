/*
 * dual_quadratic_interpolation.c
 *
 *  Created on: Dec 6, 2018
 *      Author: winona.wei
 */

#include "../include/dual_quadratic_interpolation.h"

void dual_quadratic_interpolation(double in_x[],int xn,double in_y[],int yn,double z[6][6],double h,double t,double x,double y,double *out_p)
{
	double x_i = 0;
	double y_j = 0;

	double x_l[3] = {0};
	double y_l[3] = {0};

	int k = 0;
	int r = 0;

	for(int i = 0; i < xn-1; i++)
	{
		x_i = in_x[i];
		for(int j = 0; j < yn-1; j++)
		{
			//判断x,y满足的范围
			y_j = in_y[j];

			if(1 == i && x <= (x_i + h/2) )
			{
				//用于插值的三个点下标0，1，2
				k = i - 1;
			}
			if (i >= 2 && i <= xn - 3 && x > (x_i - h/2) && x <= (x_i + h/2))
			{
				//用于插值的三个点下标i-1，i，i+1
				k = i - 1;
			}

			if(xn-2 == i && x > (x_i - h/2))
			{
				//用于插值的三个点下标n-2=x_n-3，n-1=x_n-2，n=x_n-1
				k = i - 1;
			}

			if(1 == j && y <= (y_j + t/2))
			{
				r = j - 1;
			}

			if(j >= 2 && j <= yn - 3 && y > (y_j - t/2) && y <= (y_j + t/2))
			{
				r = j - 1;
			}
			if(yn-2 == i && y > (y_j - t/2))
			{
				r = j - 1;
			}

			x_l[0] = ((x-in_x[k+1])*(x-in_x[k+2]))/((in_x[k]-in_x[k+1])*(in_x[k]-in_x[k+2]));
			x_l[1] = ((x-in_x[k])*(x-in_x[k+2]))/((in_x[k+1]-in_x[k])*(in_x[k+1]-in_x[k+2]));
			x_l[2] = ((x-in_x[k])*(x-in_x[k+1]))/((in_x[k+2]-in_x[k])*(in_x[k+2]-in_x[k+1]));

			y_l[0] = ((y-in_y[r+1])*(y-in_y[r+2]))/((in_y[r]-in_y[r+1])*(in_y[r]-in_y[r+2]));
			y_l[1] = ((y-in_y[r])*(y-in_y[r+2]))/((in_y[r+1]-in_y[r])*(in_y[r+1]-in_y[r+2]));
			y_l[2] = ((y-in_y[r])*(y-in_y[r+1]))/((in_y[r+2]-in_y[r])*(in_y[r+2]-in_y[r+1]));

			*out_p = 0;

			for(int q = k;q < k+3;q++)
			{
				for(int p = r;p < r+3;p++)
				{
					*out_p = *out_p + x_l[q-k] * y_l[p-r] * z[q][p];
				}
			}

		}
	}

	return;
}
