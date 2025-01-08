#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <complex>
#include <algorithm>
using namespace std;

/*param:
		gauss - array with all values set to zero;
		gauss_freq - array with frequences where we generate our gaussian function
		gauss_amp - amplitude of a gauss function on a given freq
*/
vector<uint16_t> gauss_freq{31, 62, 125, 250, 500, 1000, 2000, 4000, 8000, 16000};

const uint8_t BASS = 0;
const uint8_t MID = 1;
const uint8_t TREBLE = 2;

vector<vector<int>> gauss_amp{
	{4, 3, 2, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 4, 4, 4, 4, 0, 0, 0},
	{4, 3, 2, 0, 0, 0, 0, 2, 3, 4}};

vector<double> generate_gaussian(vector<double> freq, uint8_t preset)
{
	vector<double> gauss(freq.size());
	for (int i = 0; i < freq.size(); i++)
	{
		gauss[i] = 0;
	}

	for (int i = 0; i < gauss_freq.size(); i++)
	{
		for (int j = 0; j < freq.size(); j++)
			gauss[j] += gauss_amp[preset][i] * exp(-(pow(freq[j] - gauss_freq[i], 2) / (2 * pow(gauss_freq[i] / 3, 2))));
	}

	return gauss;
}

/*
Peak normalization of a fft Array.

param:
		fftArr - the vector that needs to be normalzied.
		returns vector that is normalized

*/

vector<complex<double>> normalizePeak(vector<complex<double>> fftArr, uint8_t preset)
{
	vector<complex<double>> out(fftArr.size());

	int normalizeValue = *max_element(gauss_amp[preset].begin(), gauss_amp[preset].end()) + 1;

	for (int i = 0; i < fftArr.size(); i++)
	{
		out[i].real(fftArr[i].real() / normalizeValue);
	}
	return out;
}

/*
Power normalization:

param:
		modifiedFftArr - the fft vector that we modified by the gauss array
		fftArr - the unmodofied vector
		returns vector that is normalized

*/

vector<complex<double>> normalizePower(vector<complex<double>> modifiedFftArr, vector<complex<double>> fftArr)
{
	double power1 = 0;
	double power2 = 0;
	double multiplier;

	// calculate the total power of the unmodified signal
	for (int i = 0; i <= fftArr.size(); i++)
	{
		power1 += abs(fftArr[i].real());
	}
	// calculate the total power of the modified signal
	for (int i = 0; i <= fftArr.size(); i++)
	{
		power2 += abs(modifiedFftArr[i].real());
	}

	// the square of the ratios of the 2 powers gives me a multiplier which we will use to normalize our signal
	multiplier = pow(power1 / power2, 2);
	vector<complex<double>> out(fftArr.size());

	// multiply every element of the array
	for (int i = 0; i < fftArr.size(); i++)
	{
		out[i] = modifiedFftArr[i] * multiplier;
	}

	return out;
}

// Finding nearest number power of 2
int resize(int samples)
{
	int nth_degree[21] = {1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824};

	for (int i = 0; i < 21; i++)
	{
		if (samples < nth_degree[i])
		{
			return nth_degree[i];
			break;
		}
	}
	return samples;
}
