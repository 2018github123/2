/*
 * mymain.c
 *
 *  Created on: Dec 6, 2018
 *      Author: winona.wei
 */

#include "../include/dual_quadratic_interpolation.h"
#include "../include/newton.h"
#include "../include/LS.h"
#include "../include/common.h"

int main() {

	double vec_x[4] = { 0 };
	double t[6] = { 0, 0.2, 0.4, 0.6, 0.8, 1.0 };
	double u[6] = { 0, 0.4, 0.8, 1.2, 1.6, 2.0 };

	double x[X_num] = { 0 };
	double y[Y_num] = { 0 };
	double z_xy[X_num][Y_num];
	double p = 0;
	double sigma = 0;
	double z[6][6]={{-0.5,-0.34,0.14,0.94,2.06,3.5},
					{-0.42,-0.5,-0.26,0.3,1.18,2.38},
					{-0.18,-0.5,-0.5,-0.18,0.46,1.42},
					{0.22,-0.34,-0.58,-0.5,-0.1,0.62},
					{0.78,-0.02,-0.5,-0.66,-0.5,-0.02},
					{1.5,0.46,-0.26,-0.66,-0.74,-0.5}};
	int k = 0;
	double f[X_num][Y_num];
	memset(f,0,X_num*Y_num*sizeof(double));
	memset(z_xy, 0, X_num * Y_num * sizeof(double));

	for (int j = 0; j < Y_num; j++) {
		y[j] = 0.5 + 0.05 * j;
	}
//------------------------------------
	printf("初始化向量x0=[1,1,1,1]\n");
	for (int i = 0; i < X_num; i++) {
		x[i] = 0.08 * i;
		for (int j = 0; j < Y_num; j++) {
			newton(x[i], y[j], vec_x);
			printf("(x,y)=(%f,%f)，(t,u)=(%f,%f)\t\n",x,y,vec_x[0],vec_x[1]);
			//求出x,y ->t,u->z
		//	dual_quadratic_interpolation(t, 6, u, 6, z, 0.2, 0.4, vec_x[0],
		//			vec_x[1], &z_xy[i][j]);
		//	printf("(x,y)=(%f,%f)，f(x,y)=%1.11e\t\n",x[i],y[j],z_xy[i][j]);
			//	printf("%1.11e,",x[i],y[j],z_xy[i][j]);
		}
		//printf("\n");
	}
	/*
	//利用求出f(x,y)的值求出插值函数p
	for (k=0; k<10; k++)
	{
		sigma = 0;
		for (int i = 0; i < X_num; i++)
		{
			for (int j = 0; j < Y_num; j++)
			{
				p = 0;
				LS(0.08 * i, 0.5 + 0.05 * j, z_xy, k, &p);
				//printf("p[%d,%d] is %1.11e\n", i, j, p);
				sigma += pow(p - z_xy[i][j], 2);
			}
		}
		printf("k = %d\t,sigma = % 1.11e\t\n",k,sigma);
		if(sigma <= 1.0e-7)
			{
				printf("达到精度要求的k, sigma:\n");
				printf("k = %d\n",k);
				printf("sigma = %1.11e\n",sigma);
				break;
			}
	}
	*/
//----------------------------------------
	/*
	for (int i=1;i<=8;i++)
	{
			for (int j=1;j<=5;j++)
			{
				memset(vec_x,0,4*sizeof(double));
				newton(0.1*i, 0.5+0.2*j, vec_x);
				dual_quadratic_interpolation(t, 6, u, 6, z, 0.2, 0.4, vec_x[0],
						vec_x[1], &f[i-1][j-1]);
				printf("(x,y)=(%f,%f)，f(x,y)=%1.11e\t\n",0.1*i,0.5+0.2*j,f[i-1][j-1]);
			}
		}

		for (k=5; k==5; k++)
		{
			sigma = 0;
			for (int i = 1; i <= 8; i++)
			{
				for (int j = 1; j <= 5; j++)
				{
					p = 0;
					LS(0.1*i, 0.5+0.2*j, z_xy, k, &p);
					printf("(x,y)=(%f,%f), f(x,y)=%1.11e\t, p(x,y)=%1.11e\t\n",0.1 * i, 0.5 + 0.2 * j,f[i-1][j-1],p);
				}
			}
		}
*/
	return 0;
}

