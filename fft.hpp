//	https://github.com/tkaczenko/WavReader/blob/master/WavReader/WavReader.cpp
//  https://cplusplus.com/forum/general/265589/

#define _CRT_SECURE_NO_DEPRECATE
#include <cstdint>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

const double pi = 3.14159265358979323846;
using Complex = complex<double>;

// magic
vector<Complex> fft(vector<Complex> in)
{

	vector<Complex> x(in);

	// DFT
	unsigned int N = x.size(), k = N, n;
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

	return x;
}

// inverse fft (in-place)
// magic no2
vector<Complex> ifft(vector<Complex> in)
{
	vector<Complex> x1(in);
	size_t size = in.size();
	for (int i = 0; i < size; i++)
	{
		x1[i] = std::conj(x1[i]);
		// cout << x1[i] << endl;
	}

	// forward fft
	x1 = fft(x1);

	// conjugate the complex numbers again
	for (int i = 0; i < size; i++)
		x1[i] = std::conj(x1[i]);

	for (int i = 0; i < size; i++)
		x1[i] /= size;

	return x1;
}

/*
returns a vector with the frequencies

param: 
		n : length of the soundfile
		sampleRate: self Explanatory
*/
vector<double> getFreq(int n, int sampleRate)
{
	vector<double> freq(n);
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
	return freq;
}
