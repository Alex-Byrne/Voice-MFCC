
#include <stdint.h>
#include <math.h>
#include "DSPF_sp_fftSPxSP_cn.h"
#include <stdio.h>

#include <string.h>
#include "MFCC.h"

// Global Definitions and Variables
#define PI 3.14159265358979323

//#define d_100 16000

#define Nyquist 8000.0
#define ALPHA ((float)0.54)
#define BETA ((float)0.46)


#pragma DATA_ALIGN(x_sp,8);
float   x_sp[2 * N];
#pragma DATA_ALIGN(y_sp,8);
float   y_sp[2 * N];
#pragma DATA_ALIGN(w_sp,8);
float   w_sp[2 * N];
#pragma DATA_ALIGN(w_sp,8);
float lambda[2 * N];
#pragma DATA_ALIGN(w_sp,8);
float hamming[N];   

/* Align the tables that we have to use */

// The DATA_ALIGN pragma aligns the symbol in C, or the next symbol declared in C++, to an alignment boundary.
// The alignment boundary is the maximum of the symbol's default alignment value or the value of the constant in bytes.
// The constant must be a power of 2. The maximum alignment is 32768.
// The DATA_ALIGN pragma cannot be used to reduce an object's natural alignment.

//The following code will locate mybyte at an even address.
//#pragma DATA_ALIGN(mybyte, 2)
//char mybyte;

//The following code will locate mybuffer at an address that is evenly divisible by 1024.
//#pragma DATA_ALIGN(mybuffer, 1024)
//char mybuffer[256];


// brev routine called by FFT routine
unsigned char brev[64] = {
	0x0, 0x20, 0x10, 0x30, 0x8, 0x28, 0x18, 0x38,
	0x4, 0x24, 0x14, 0x34, 0xc, 0x2c, 0x1c, 0x3c,
	0x2, 0x22, 0x12, 0x32, 0xa, 0x2a, 0x1a, 0x3a,
	0x6, 0x26, 0x16, 0x36, 0xe, 0x2e, 0x1e, 0x3e,
	0x1, 0x21, 0x11, 0x31, 0x9, 0x29, 0x19, 0x39,
	0x5, 0x25, 0x15, 0x35, 0xd, 0x2d, 0x1d, 0x3d,
	0x3, 0x23, 0x13, 0x33, 0xb, 0x2b, 0x1b, 0x3b,
	0x7, 0x27, 0x17, 0x37, 0xf, 0x2f, 0x1f, 0x3f
};



// The seperateRealImg function separates the real and imaginary data
// of the FFT output. This is needed so that the data can be plotted
// using the CCS graph feature

// Function for generating sequence of twiddle factors
void gen_twiddle_fft_sp()
{
	int i, j, k;
	double x_t, y_t, theta1, theta2, theta3;

	for (j = 1, k = 0; j <= N >> 2; j = j << 2){
		for (i = 0; i < N >> 2; i += j){
			theta1 = 2 * PI * i / N;
			x_t = cos(theta1);
			y_t = sin(theta1);
			w_sp[k] = (float)x_t;
			w_sp[k + 1] = (float)y_t;

			theta2 = 4 * PI * i / N;
			x_t = cos(theta2);
			y_t = sin(theta2);
			w_sp[k + 2] = (float)x_t;
			w_sp[k + 3] = (float)y_t;

			theta3 = 6 * PI * i / N;
			x_t = cos(theta3);
			y_t = sin(theta3);
			w_sp[k + 4] = (float)x_t;
			w_sp[k + 5] = (float)y_t;
			k += 6;
		}
	}
}

void gen_hamming(){
	int n;

	for (n = 0; n < N; n++)
		hamming[n] = ALPHA - BETA * cos(2 * PI * n / (double)(N - 1));

}

void gen_lambda(){
	int n;
	int fmel[M + 2];
	float Phi_K[M + 2];

	for (n = 0; n < M + 2; n++) 
		fmel[n] = 344 + (2840 - 344) * n / 27;

	for (n = 0; n < M + 2; n++)
		Phi_K[n] = (pow(10.0, (fmel[n] / 2595.0)) - 1.0) * 700.0;

	for (n = 0; n < M + 2; n++)
		lambda[n] = (double)(N / 2 + 1) * Phi_K[n] / Nyquist;
}

float input_adj[N];
float sigma[M];
float mag[N];
float pInterest[K];
float filtered[N];
float Hk[M][K];

void gen_MFCC(int* input, int nInputs, float* mfcc)
{
	int n;

	//Zero padding (probably not best method, but is the simplest) 
	for (n = 0; n < N; n++)
		if (n > nInputs)
			input_adj[n] = 0;
		else
			input_adj[n] = input[n];

	//Hamming window (Step 1)
	for (n = 0; n < N; n++)
		filtered[n] = hamming[n] * input_adj[n];

	for (n = 0; n < N; n++){
		x_sp[2 * n] = filtered[n];
		x_sp[2 * n + 1] = 0;
	}
	//FFt (Step 2)
	DSPF_sp_fftSPxSP_cn(N, x_sp, w_sp, y_sp, brev, 4, 0, N);

	//|FFT|^2
	for (n = 0; n < N; n++)
		mag[n] = (pow(y_sp[n], 2) + pow(y_sp[n + 1], 2));

	//Find first K points
	for (n = 0; n < K; n++)
		pInterest[n] = mag[n];

	//Filter bank (Step 3)

	int k;
	for (n = 1; n < M + 2; n++) {
		for (k = 0; k < K ; k++) {
			if (k < lambda[n - 1])
				Hk[n - 1][k] = 0;
			else if (k >= lambda[n - 1] && k <= lambda[n])
				Hk[n - 1][k] = (float)(k - lambda[n - 1]) / (float)(lambda[n] - lambda[n - 1]);
			else if (k >= lambda[n] && k <= lambda[n + 1])
				Hk[n - 1][k] = (float)(lambda[n + 1] - k) / (float)(lambda[n + 1] - lambda[n]);
			else if (k > lambda[n + 1])
				Hk[n - 1][k] = 0;
			else
				printf("Error");

		}
	}

	float sum = 0;

	for (n = 0; n < M; n++) {
		for (k = 0; k < K; k++) {
			sum += (float)pInterest[k] * (float)Hk[n][k];
		}
		sigma[n] = log10f(sum);
		sum = 0;
	}

	sum = 0;
	//DCT
	int j;
	for (n = 0; n < 13; n++) {
		for (j = 0; j < M; j++) {
			sum += sigma[j] * (float)cos(((float)j - 0.5) * (float)n * PI / (float)M);
		}
		mfcc[n] = sum;
		sum = 0;
	}
}

void init_MFCC()
{
	gen_twiddle_fft_sp();
	gen_hamming();
	gen_lambda();
}

