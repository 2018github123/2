/*
 * newton.c
 *
 *  Created on: Dec 6, 2018
 *      Author: winona.wei
 */
#include "../include/newton.h"
#include "../include/gauss.h"
#include <math.h>
#include <stdio.h>

#define  M  500

void newton(double x,double y,double out_x[4])
{
	double F[4] = {0};
	double dF[4][4];

	double ksi[4] = {0};
	double detx[4] = {0};

	double max_ksi = 0;
	double max_detx = 0;

	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			dF[i][j] = 1;
		}
	}

	dF[2][0] = 0.5;
	dF[3][1] = 0.5;

	//在定义域范围内初始化ksi
	// ksi[4]
	for(int i= 0;i<4;i++)
	{
		ksi[i]=1000.0;
		out_x[i] = 0;
	}

//循环M次
	for(int k=0;k<M;k++)
	{
		//计算F dF
		F[0] = -1.0 * (0.5*cos(ksi[0])+ksi[1]+ksi[2]+ksi[3]-x-2.67);
		F[1] = -1.0 * (ksi[0]+0.5*sin(ksi[1])+ksi[2]+ksi[3]-y-1.07);
		F[2] = -1.0 * (0.5*ksi[0]+ksi[1]+cos(ksi[2])+ksi[3]-x-3.74);
		F[3] = -1.0 * (ksi[0]+0.5*ksi[1]+ksi[2]+sin(ksi[3])-y-0.79);

		dF[0][0]= -0.5 * sin(ksi[0]);
		dF[1][1]= 0.5 * cos(ksi[1]);
		dF[2][2]= -1.0 * sin(ksi[2]);
		dF[3][3]= cos(ksi[3]);
		//解线性方程
		for(int i=0;i<4;i++)
		{
			detx[i] = 0;
		}
		gauss(dF,F,detx);
		//
		max_ksi = fabs(ksi[0]);
		max_detx = fabs(detx[0]);

		for(int i=1;i<4;i++)
		{
			if(fabs(ksi[i]) > max_ksi)
			{
				max_ksi = fabs(ksi[i]);
			}
			if(fabs(detx[i]) > max_detx)
			{
				max_detx = fabs(detx[i]);
			}
		}
		//精度判断
		if((max_detx/max_ksi) > 1.0e-12)
		{
			ksi[0] = ksi[0] + detx[0];
			ksi[1] = ksi[1] + detx[1];
			ksi[2] = ksi[2] + detx[2];
			ksi[3] = ksi[3] + detx[3];

			if(k == M)
			{
				printf("非线性方程组迭代不成功！\n");
			}
		}
		else
		{
			for(int i= 0;i<4;i++)
			{
				out_x[i] = ksi[i];
			}
			break;
		}
	}

	return;
}
