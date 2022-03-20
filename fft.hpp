//	https://github.com/tkaczenko/WavReader/blob/master/WavReader/WavReader.cpp
//  https://cplusplus.com/forum/general/265589/

#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <algorithm>		// For std::transform
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
#include "string.h"
#include "memory.h"

using namespace std;


const double pi = 3.14159265358979323846;
using Complex = complex<double>;

// Fast Fourier Transform
void FFT(Complex f[], Complex ftilde[], int log2N)
{
	int N = 1 << log2N;
	cout << "N: " <<N<<endl;
	cout << "LOG2N: " << log2N<<endl;
	int iii;
	// Reorder
	for (int i = 0; i < N; i++)
	{
		int ii = 0, x = i;
		for (int j = 0; j < log2N; j++)
		{
			ii <<= 1;
			ii |= (x & 1);
			x >>= 1;
		}
		ftilde[ii] = f[i];
		iii = ii;
	}

	cout  << "ii: " << iii << endl;

	for (int s = 1; s <= log2N; s++)
	{
		int m = 1 << s;
		int m2 = m / 2;
		Complex w = 1.0;
		Complex wm = polar(1.0, -pi / m2);
		for (int j = 0; j < m2; j++)
		{
			for (int k = j; k < N; k += m)
			{
				Complex t = w * ftilde[k + m2];
				Complex u = ftilde[k];
				ftilde[k] = u + t;
				ftilde[k + m2] = u - t;
			}
			w *= wm;
		}
		
	}

}

// Inverse Fast Fourier Transform
void iFFT(Complex ftilde[], Complex f[], int log2N)
{
	int N = 1 << log2N;

	for (int m = 0; m < N; m++) ftilde[m] = conj(ftilde[m]);      // Apply conjugate (reversed below)

	FFT(ftilde, f, log2N);

	double factor = 1.0 / N;
	for (int m = 0; m < N; m++) f[m] = conj(f[m]) * factor;

	for (int m = 0; m < N; m++) ftilde[m] = conj(ftilde[m]);      // Only necessary to reinstate ftilde
}

void fft(Complex x[],int size)
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
 