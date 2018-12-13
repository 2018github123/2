/*
 * QR.c
 *
 *  Created on: Dec 7, 2018
 *      Author: winona.wei
 *
 */

#include "../include/QR.h"

//input A  (k+1)*(k+1)

void QR(double A[NUM][NUM],double out_Q[NUM][NUM],int n)
{
	int m = n-1;
	double sum = 0;
	double dr = 0;
	double cr = 0;
	double hr = 0;

	double Q[NUM][NUM];

	double ur[NUM] = {0};
	double pr[NUM] = {0};
	double wr[NUM] = {0};
	//初始化
	for (int i = 0; i < NUM; i++) {
		for (int j = 0; j < NUM; j++) {
			Q[i][j] = 0;
			out_Q[i][j] = 0;
		}

	}


	for(int i = 0;i<n;i++)
	{
		Q[i][i] = 1;
	}

	// Mk QR分解算法
	//r=0,1,2,...,n-2
	for(int r = 0; r < m ; r++)
	{
		//(1)
		sum = 0;
		for(int i = r + 1; i <= m; i++)
		{
			sum = sum + fabs(A[i][r]);
		}
		if(0 == sum)
		{
			continue;
		}

		//(2)
		dr = 0;
		for(int i = r;i <= m; i++)
		{
			dr = dr + A[i][r] *  A[i][r];
		}
		dr = sqrt(dr);

		if(A[r][r] > 0 )
		{
			cr = -1.0 * dr;
		}
		if(A[r][r] < 0 ||  A[r][r] == 0)
		{
			cr = dr;
		}
		hr = cr * cr - cr * A[r][r];
		//(3)

		for (int i = 0; i < NUM; i++) {
				ur[i] = 0;
				pr[i] = 0;
				wr[i] = 0;
		}

		ur[r] = A[r][r]-cr;
		for (int i = r+1; i <= m; i++) {
			ur[i] = A[i][r];
		}
		//(4)

		for(int i=0;i<=m;i++)
		{
			for(int j=0;j<=m;j++)
			{
				wr[i] += Q[i][j] * ur[j];
				pr[i] += A[j][i] * ur[j];
			}
			pr[i] /= hr;
		}


		for(int i=0;i<=m;i++)
		{
			for(int j=0;j<=m;j++)
			{
				Q[i][j] = Q[i][j] - wr[i]*ur[j]/hr;
			}
		}

		//A(r+1)
		for(int i = 0; i <= m; i++)
		{
			for(int j = 0;j <= m; j++)
			{
				A[i][j] = A[i][j] - ur[i]*pr[j];//An = R
			}
		}
	}

	for(int i = 0; i <= m; i++)
	{
		for(int j = 0;j <= m; j++)
		{
			out_Q[i][j] = Q[i][j];
		}
	}

	return;
}

