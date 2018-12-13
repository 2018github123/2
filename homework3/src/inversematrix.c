#include "../include/inversematrix.h"
#include "../include/common.h"

void inversematrix(double A[NUM][NUM], double inv_A[NUM][NUM], int k) {
	double ATA[NUM][NUM];
	double Q[NUM][NUM];
	double QT[NUM][NUM];
	double inR[NUM][NUM];

	for (int i = 0; i < NUM; i++) {
		for (int j = 0; j < NUM; j++) {
			ATA[i][j] = 0;
			Q[i][j] = 0;
			QT[i][j] = 0;
			inR[i][j] = 0;
			inv_A[i][j] = 0;
		}

	}

	for (int i = 0; i < (k + 1); i++) {
		inR[i][i] = 1;
	}

	for (int i = 0; i < (k + 1); i++) {
		for (int j = 0; j < (k + 1); j++) {
			ATA[i][j] = A[i][j];
		}
	}

	//｛ATA｝ (k+1)*(k+1)的逆
	//QR分解
	QR(ATA, Q, k + 1);

	for (int j = k; j >= 0; j--) {
		for (int t = j; t <= k; t++) {
			inR[j][t] = inR[j][t] / ATA[j][j];	//对角线ATA化简为1
		}
		for (int i = j - 1; i >= 0; i--) {
			for (int t = j; t <= k; t++) {
				inR[i][t] -= ATA[i][j] * inR[j][t];
			}
		}
	}

	for (int i = 0; i < (k + 1); i++) {
		for (int j = 0; j < (k + 1); j++) {
			QT[i][j] = Q[j][i];
		}
	}
	double sum = 0;
	for (int i = 0; i < (k + 1); i++) {
		for (int j = 0; j < (k + 1); j++) {
			sum = 0;
			for (int t = 0; t < (k + 1); t++) {
				sum += inR[i][t] * QT[t][j];
			}
			inv_A[i][j] = sum;
		}
	}


	return;
}

