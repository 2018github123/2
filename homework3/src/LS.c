/*
 * LS.c
 *
 *  Created on: Dec 7, 2018
 *      Author: yanpengfei
 */

#include "../include/LS.h"
#include "../include/inversematrix.h"

void LS(double x, double y, double u[X_num][Y_num], int k, double *p) {
	double B[NUM][NUM];
	double G[NUM][NUM];
	double B_T[NUM][NUM];
	double G_T[NUM][NUM];
	double BTB[NUM][NUM];
	double GTG[NUM][NUM];

	double inv_BB[NUM][NUM];
	double inv_GG[NUM][NUM];

	double tmp1[NUM][NUM];
	double tmp2[NUM][NUM];
	double tmp3[NUM][NUM];

	double U[X_num][Y_num];
	double c[X_num][Y_num];
	//printf("U :\n");
	for (int i = 0; i < X_num; i++) {
		for (int j = 0; j < Y_num; j++) {
			U[i][j] = u[i][j];
			c[i][j] = 0;
	//		printf("% 1.11e\t",U[i][j]);
		}
	//	printf("\n");
	}

	for (int i = 0; i < NUM; i++) {
			for (int j = 0; j < NUM; j++) {
				B[i][j] = 0;
				G[i][j] = 0;
				B_T[i][j] = 0;
				G_T[i][j] = 0;
				BTB[i][j]=0;
				GTG[i][j]=0;
				inv_BB[i][j] = 0;
				inv_GG[i][j] = 0;
				tmp1[i][j] = 0;
				tmp2[i][j] = 0;
				tmp3[i][j] = 0;
			}
		}
	//printf("B:\n");
	for (int i = 0; i < X_num; i++) {
		for (int r = 0; r < (k + 1); r++) {
			B[i][r] = pow(0.08 * i, r);
	//		printf("% 1.11e\t",B[i][r]);
		}
	//	printf("\n");
	}

	//printf("B_T:\n");
	for (int i = 0; i < (k + 1); i++) {
		for (int j = 0; j < X_num; j++) {
			B_T[i][j] = B[j][i];
	//		printf("% 1.11e\t",B_T[i][j]);
		}
	//	printf("\n");
	}

	//printf("BTB:\n");

	for (int i = 0; i < (k + 1); i++)
	{
		for (int j = 0; j < (k + 1); j++)
		{
			for (int t = 0; t < X_num; t++)
			{
				//A的转置＊A
				BTB[i][j] += B_T[i][t] * B[t][j];

			}
	//		printf("% 1.11e\t",BTB[i][j]);
		}
	//	printf("\n");
	}

//	printf("\n");
//	printf("G :\n");
	for (int i = 0; i < Y_num; i++)
	{
		for (int s = 0; s < (k + 1); s++)
		{
			G[i][s] = pow(0.5 + 0.05 * i, s);
	//		printf("% 1.11e\t",G[i][s]);
		}
	//	printf("\n");
	}

	//printf(" G_T\n");
	for (int i = 0; i < (k + 1); i++)
	{
		for (int j = 0; j < Y_num; j++)
		{
			G_T[i][j] = G[j][i];
	//		printf("% 1.11e\t",G_T[i][j]);

		}
	//	printf("\n");
	}
	//printf(" GTG\n");
	for (int i = 0; i < (k + 1); i++)
	{
		for (int j = 0; j < (k + 1); j++)
		{
			for (int t = 0; t < Y_num; t++) {
				//A的转置＊A
				GTG[i][j] += G_T[i][t] * G[t][j];
			}
	//		printf("% 1.11e\t",GTG[i][j]);
		}
	//	printf("\n");
	}
	inversematrix(BTB, inv_BB, k);
	inversematrix(GTG, inv_GG, k);
	//printf("inv_BB:\n");




	//c = inv_BB*B_T*u*G*inv_GG;

	//tmp1 = inv_BB*B_T
	//printf("tmp1:\n");
	for (int i = 0; i < (k + 1); i++)
	{
		for (int j = 0; j < X_num; j++)
		{
			for (int t = 0; t < (k + 1); t++)
			{
				tmp1[i][j] += inv_BB[i][t] * B_T[t][j];
			}
		//	printf("%1.11e\t",tmp1[i][j]);
		}
		//printf("\n");
	}
//	printf("tmp2:\n");
	for (int i = 0; i < (k + 1); i++)
	{
		for (int j = 0; j < Y_num; j++)
		{
			for (int t = 0; t < X_num; t++)
			{
				tmp2[i][j] += tmp1[i][t] * U[t][j];
			}
			//printf("%1.11e\t",tmp2[i][j]);
		}
	//	printf("\n");
	}
	//printf("tmp3:\n");
	for (int i = 0; i < (k + 1); i++)
	{
		for (int j = 0; j < (k + 1); j++)
		{
			for (int t = 0; t < Y_num; t++)
			{
				tmp3[i][j] += tmp2[i][t] * G[t][j];
			}
			//printf("%1.11e\t",tmp3[i][j]);
		}
		//printf("\n");
	}
	printf("c is :\n");
	for (int i = 0; i < (k + 1); i++)
	{
		for (int j = 0; j < (k + 1); j++)
		{
			for (int t = 0; t < (k + 1); t++)
			{
				c[i][j] += tmp3[i][t] * inv_GG[t][j];
			}
			printf("% 1.11e\t",c[i][j]);
		}
		printf("\n");
	}

	//计算精度，确定最小的k
	/*printf("ls x is %f\n",x);
	printf("ls y is %f\n",y);*/
	*p = 0;
	for (int r = 0; r < k + 1; r++)
	{
		for (int s = 0; s < k + 1; s++)
		{
			(*p) += c[r][s] * pow(x, r) * pow(y, s);
			//printf("c[r][s]=%1.11e\t",c[r][s]);
		}
		//printf("\n");
	}

	//printf("ls p is %1.11e\n",*p);

	return;
}

