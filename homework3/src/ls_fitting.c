/*

 * ls_fitting.c
 *
 *  Created on: Dec 6, 2018
 *      Author: winona.wei


#include "../include/ls_fitting.h"
#include "../include/QR.h"

#define M 10 //迭代次数  = min(X_num,Y_num)


void inversematrix(double A[NUM][NUM],double inv_A[NUM][NUM],int row_num,int k);

void ls_fitting(double x[X_num],double y[Y_num],const double u[X_num][Y_num],double c[X_num][Y_num],int *out_k,double *out_sigma)
{
	double B[NUM][NUM];
	double G[NUM][NUM];
	double B_T[NUM][NUM];
	double G_T[NUM][NUM];
	double BTB[NUM][NUM];
	double GTG[NUM][NUM];

	double inv_BB[NUM][NUM];
	double inv_GG[NUM][NUM];

	double p[NUM][NUM];
	double sigma = 0;

	double tmp1[NUM][NUM];
	double tmp2[NUM][NUM];
	double tmp3[NUM][NUM];
	double U[X_num][Y_num];

	//
	//printf("u is :\n");

	for(int i = 0;i<X_num;i++)
	{

		for(int j = 0;j<Y_num;j++)
		{
			U[i][j] = u[i][j];
			//printf("% 1.11e\t",U[i][j]);
		}
		//printf("\n");
	}

	for(int k=0;k<M;k++)
	{
		//初始化为0
		for (int i = 0; i < NUM; i++) {
			for (int j = 0; j < NUM; j++) {
				B[i][j] = 0;
				G[i][j] = 0;
				inv_BB[i][j] = 0;
				inv_GG[i][j] = 0;
				B_T[i][j] = 0;
				p[i][j] = 0;
				tmp3[i][j] = 0;
				tmp3[i][j] = 0;
				tmp3[i][j] = 0;
				BTB[i][j]=0;
				GTG[i][j]=0;
			}

		}

		for(int i=0;i<X_num;i++)
		{
			for(int j=0;j<Y_num;j++)
			{
				c[i][j] = 0;
			}
		}

		//printf(" B\n");
		//固定y
		for(int i = 0;i<X_num;i++)
		{
			for(int r=0;r<(k+1);r++)
			{
				B[i][r] = pow(x[i],r);

			//	printf("B[%d][%d]=% 1.11e\t",i,r,B[i][r]);
			}
			//printf("\n");
		}
		//printf(" B_T\n");
		for(int i=0;i<(k+1);i++)
		{
			for(int j = 0;j<X_num;j++)
			{
				B_T[i][j] = B[j][i];
				//printf("B_T[%d][%d]=% 1.11e\t",r,i,B_T[r][i]);
			}
			//printf("\n");
		}
		for(int i=0;i<(k+1);i++)
		{
			for(int j=0;j<(k+1);j++)
			{
				BTB[i][j] = 0;
				for(int p=0;p<X_num;p++)
				{
					//A的转置＊A
					BTB[i][j] += B_T[i][p]*B[p][j];
				}
			}
		}


		for(int i = 0;i<Y_num;i++)
		{
			for(int s=0;s<(k+1);s++)
			{
				G[i][s] = pow(0.5+0.05*i,s);
			}
		}

		//printf(" B_T\n");
		for (int i = 0; i < (k + 1); i++) {
			for (int j = 0; j < Y_num; j++) {
				G_T[i][j] = G[j][i];
				//printf("B_T[%d][%d]=% 1.11e\t",r,i,B_T[r][i]);
			}
			//printf("\n");
		}

		for (int i = 0; i < (k + 1); i++) {
			for (int j = 0; j < (k + 1); j++) {
				GTG[i][j] = 0;
				for (int p = 0; p < Y_num; p++) {
					//A的转置＊A
					GTG[i][j] += G_T[i][p] * G[p][j];
				}
			}
		}

		inversematrix(BTB,inv_BB,X_num,k);
		inversematrix(GTG,inv_GG,Y_num,k);


		//c = inv_BB*B_T*u*G*inv_GG;

		//tmp1 = inv_BB*B_T
		for (int i = 0; i < (k+1); i++) {
			for (int j = 0; j < X_num; j++) {
				for(int t = 0; t < (k+1); t++)
				{
					tmp1[i][j] += inv_BB[i][t] * B_T[t][j];
				}
			}
		}

		for (int i = 0; i < (k+1); i++) {
			for (int j = 0; j < Y_num; j++) {
				for(int t = 0; t < X_num; t++)
				{
					tmp2[i][j] += tmp1[i][t] * U[t][j];
				}
			}
		}


		for (int i = 0; i < (k + 1); i++) {
			for (int j = 0; j < (k + 1); j++) {
				for (int t = 0; t < Y_num; t++) {
					tmp3[i][j] += tmp2[i][t] * G[t][j];
				}
			}
		}
		//printf("c is :\n");
		for (int i = 0; i < (k + 1); i++) {
			for (int j = 0; j < (k + 1); j++) {
				for (int t = 0; t < (k + 1); t++) {
					c[i][j] += tmp3[i][t] * inv_GG[t][j];
				}
				//printf("% 1.11e\t",c[i][j]);
			}
			//printf("\n");
		}

		//计算精度，确定最小的k
		sigma = 0;
		//printf("U is :\n");
		for(int i=0;i<X_num;i++)
		{
			for(int j=0;j<Y_num;j++)
			{
				p[i][j] = 0;
				for(int r=0;r<k+1;r++)
				{
					for(int s=0;s<k+1;s++)
					{
						p[i][j] = p[i][j]+c[r][s]*pow(x[i],r)*pow(y[j],s);

					}
				}

				sigma += pow(U[i][j] - p[i][j],2);
				printf("sigma = %1.11e\t",sigma);
			}
			//printf("\n");
		}

		//printf("sigma = %1.11e\n",sigma);

		if( sigma <= 1.0e-7)
		{
			printf(" 第k=%d次循环计算\n",k);
			*out_k = k;
			*out_sigma = sigma;
			printf(" 第*out_k=%d次循环计算\n",*out_k);
			break;
		}
		else if(k == M)
		{
			printf(" failed\n");
		}
	}

	return;
}

void inversematrix(double A[NUM][NUM],double inv_A[NUM][NUM],int row_num,int k)
{
	double ATA[NUM][NUM];
	double Q[NUM][NUM];
	double QT[NUM][NUM];
	double inR[NUM][NUM];

	for (int i = 0; i < NUM; i++)
	{
		for (int j = 0; j < NUM; j++)
		{
			ATA[i][j] = 0;
			Q[i][j] = 0;
			QT[i][j] = 0;
			inR[i][j] = 0;
			inv_A[i][j] = 0;
		}

	}

	for(int i=0;i<(k+1);i++)
	{
		inR[i][i] = 1;
	}



	for (int i = 0; i < (k+1); i++) {
		for (int j = 0; j < (k+1); j++) {
			ATA[i][j] = A[i][j];
		}
	}

	printf("ATA :\n");
	for(int i=0;i<=k;i++)
	{
		for(int j=0;j<=k;j++)
		{
			printf("%1.11e\t",ATA[i][j]);
		}
		printf("\n");
	}

	printf("inR0 :\n");
	for(int i=0;i<=k;i++)
	{
		for(int j=0;j<=k;j++)
		{
			printf("%1.11e\t",inR[i][j]);
		}
		printf("\n");
	}


	//｛ATA｝ (k+1)*(k+1)的逆
	//QR分解
	QR(ATA,Q,k+1);

	for(int j=k;j>=0;j--)
	{
		for(int t=j;t<=k;t++)
		{
			inR[j][t] = inR[j][t]/ATA[j][j];//对角线ATA化简为1
		}
		for(int i=j-1;i>=0;i--)
		{
			for(int t=j;t<=k;t++)
			{
				inR[i][t] -= ATA[i][j]*inR[j][t];
			}
		}
	}
	printf("inR :\n");
	for(int i=0;i<=k;i++)
	{
		for(int j=0;j<=k;j++)
		{
			printf("%1.11e\t",inR[i][j]);
		}
		printf("\n");
	}

	for (int i = 0; i < (k+1); i++)
	{
			for (int j = 0; j <(k+1); j++)
			{
				QT[i][j] = Q[j][i];
			}
		}
	double sum = 0;
	for (int i = 0; i < (k+1); i++) {
		for (int j = 0; j < (k+1); j++) {
			sum = 0;
			for(int t = 0; t < (k+1); t++)
			{
				sum += inR[i][t] * QT[t][j];
			}
			inv_A[i][j]  = sum;
		}
	}

	return;
}



*/
