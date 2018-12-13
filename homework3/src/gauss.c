/*
 * gauss.c
 *
 *  Created on: Dec 6, 2018
 *      Author: winona.wei
 */
#include "../include/gauss.h"
#include <math.h>
#include <stdio.h>

void gauss(double AX[4][4],double bx[4],double out_x[4])
{
	//int ik = 0;
	double a_ik = 0;
	int max_i = 0;
	int max_j = 0;
	double tmp_a  = 0;
	double tmp_b = 0;
	double m_ik = 0;
	double tmp_sum = 0;
	double A[4][4];
	double b[4] = {0};

	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			A[i][j] = AX[i][j];
		}
	}

	for(int i=0;i<4;i++)
	{
		b[i]=bx[i];
		out_x[i] = 0;
	}


	for(int ik=0 ; ik < 3; ik++)
	{
		//选择行号ik=0,1,2
		a_ik = fabs(A[ik][ik]);
		max_i = ik;
		max_j = ik;
		for(int i = ik + 1; i < 4; i++)
		{
			if(fabs(A[i][ik]) > a_ik)
			{
				a_ik = fabs(A[i][ik]);
				max_i = i;
			}
		}
		//首先判断A是否非奇异
		if(a_ik == 0)
		{
			printf("Gauss failed...\n");
			return;
		}
		//
		for(int j=ik;j<4;j++)
		{
			tmp_a = A[ik][j];
			A[ik][j] =  A[max_i][j];
			A[max_i][j] =  tmp_a;
		}

		//
		tmp_b = b[ik];
		b[ik] = b[max_i];
		b[max_i] = tmp_b;
		//
		for(int i = ik + 1;i < 4;i++)
		{
			m_ik = A[i][ik] / A[ik][ik];

			for(int j = ik; j < 4;j++)
			{
				A[i][j]=A[i][j]-m_ik*A[ik][j];
			}
			b[i] = b[i] - m_ik * b[ik];
		}
	}

	//	//回代过程
	out_x[3] = b[3]/A[3][3];

	for(int k=2;k>=0;k--)
	{
		tmp_sum = 0;

		for(int j=k+1;j<4;j++)
		{
			tmp_sum += A[k][j]*out_x[j];
		}

		out_x[k] = (b[k]-tmp_sum)/A[k][k];
	}

	return;
}
