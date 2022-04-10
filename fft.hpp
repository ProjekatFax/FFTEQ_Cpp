//	https://github.com/tkaczenko/WavReader/blob/master/WavReader/WavReader.cpp
//  https://cplusplus.com/forum/general/265589/

#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <algorithm> // For std::transform
#include <string>
#include <istream>
#include <fstream>
#include <cstdint>
#include <stdio.h>
#include <vector>
#include <chrono>
#include <iomanip>
#include <complex>
#include <cstdio>
#include <cmath>
#include <math.h>
#include "string.h"
#include "memory.h"

using namespace std;

const double pi = 3.14159265358979323846;
using Complex = complex<double>;

void fft(Complex x[], int size)
{
	// DFT
	unsigned int N = size, k = N, n;
	double thetaT = 3.14159265358979323846264338328L / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				Complex t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
}

// inverse fft (in-place)
void ifft(Complex x[], Complex output[], int n)
{
	vector<Complex> x1;

	for (int i = 0; i < n; i++)
		x1.push_back(x[i]);

	for (int i = 0; i < n; i++){
		x1[i] = std::conj(x1[i]);
		//cout << x1[i] << endl;
	}

	// forward fft
	fft(x1.data(), n);

	// conjugate the complex numbers again
	for (int i = 0; i < n; i++)
		x1[i] = std::conj(x1[i]);

	for (int i = 0; i < n; i++)
		x1[i] /= n;

	output = x1.data();
}

void getFreq(double freq[], int n, int sampleRate)
{
	double d = 1.0 / sampleRate;
	double val = 1.0 / (n * d);
	int N = (n - 1) / 2 + 1;

	cout << "d: " << d << "   val: " << val << "   N = " << N << endl;

	for (int i = 0; i <= N; i++)
	{
		freq[i] = i * val;
	}

	int backWards = n / 2;
	for (int i = N + 1; i <= n; i++)
	{
		freq[i] = -backWards * val;
		backWards--;
	}
}

/*
void idft(Complex x[], Complex output[], int n)
{

	for (int i = 0; i < n; i++)
	{
		for (k = 0; k < n; k++)
		{
			int theta = (2 * 3.141592 * k * n) / N;
			x[n] = x[n] + Xr[k] * cos(theta) + Xi[k] * sin(theta);
		}
		x[n] = x[n] / N;
	}
}
*/
