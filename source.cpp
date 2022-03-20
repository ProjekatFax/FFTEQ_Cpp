//	https://github.com/tkaczenko/WavReader/blob/master/WavReader/WavReader.cpp
//  https://cplusplus.com/forum/general/265589/
//  https://github.com/adamstark/AudioFile/blob/master/README.md

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


#include "AudioFile.h"								// Include Stark
AudioFile<double> wav;								// Stark
//wav.load("C:/Users/M�sz�ros Arnold/Documents/Visual Studio 2015/Projects/Wav_with_Stark/Wav_with_Stark.wav");

using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::fstream;
using std::string;
using std::getline;

const double pi = 3.14159265358979323846;
using Complex = complex<double>;

struct  wav_header_t
{
	char chunkID[4]; 				// "RIFF" = 0x46464952
	unsigned long chunkSize; 		// 28 [+ sizeof(wExtraFormatBytes) + wExtraFormatBytes] + sum(sizeof(chunk.id) + sizeof(chunk.size) + chunk.size)
	char format[4]; 				// "WAVE" = 0x45564157
	char subchunk1ID[4]; 			// "fmt " = 0x20746D66

	unsigned long subchunk1Size; 	// 16 [+ sizeof(wExtraFormatBytes) + wExtraFormatBytes]
	unsigned short audioFormat;		// Audio format 1=PCM,6=mulaw,7=alaw,     257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
	unsigned short numChannels;		// Number of channels 1=Mono 2=Sterio
	unsigned long sampleRate;		// Sampling Frequency( in Hz )
	unsigned long byteRate;			// 
	unsigned short blockAlign;		// 2=16-bit mono, 4=16-bit stereo
	unsigned short bitsPerSample;	//
};

struct chunk_t
{
	char ID[4];						// "data" = 0x61746164
	unsigned long size;				// Chunk data bytes
};

// Function prototypes																			
void WavOpen(const char* fileName, const char* writeName, int mode, vector<double> &V, vector<complex<double>> &iFFT);			// Unified function for read and write
int resize(int samples);
void print(const char * prompt, Complex A[], int log2N);
void FFT(Complex f[], Complex ftilde[], int log2N);
void iFFT(Complex ftilde[], Complex f[], int log2N);

int main(int argc, char* argv[])
{
	vector<double> amplitudes;																											// Real points
	vector<complex<double>> input;																										// Just for initialisation WavOpen - Read function

	wav.load("input.wav");
	//wav.load("Wav_with_Stark.wav");

																																	//	WavOpen("C:/Users/M�sz�ros Arnold/Documents/visual studio 2015/Projects/Project4/input.wav", "output name", 0, amplitudes, input);	// Gives back data from wav
	WavOpen("input.wav", "output name", 0, amplitudes, input);	// Gives back data from wav

	size_t size = amplitudes.size();																									// Number equals to samples
	const int new_size = resize(size);																									// Nearest number power of 2
	int log2N = log2(new_size);																											// log2 base of nearest number
	cout << "Number of Vector elements: " << size << ", numer of new elements: " << new_size << ", with log = " << log2N << '\n';		// Checking values
	cout << '\n';

	cout << "Stark algorithm data: " << '\n';														// Checking values Stark values
	cout << '\n';
	int sampleRate = wav.getSampleRate();
	int bitDepth = wav.getBitDepth();
	int numSamples = wav.getNumSamplesPerChannel();
	double lengthInSeconds = wav.getLengthInSeconds();
	int numChannels = wav.getNumChannels();
	bool isMono = wav.isMono();
	bool isStereo = wav.isStereo();
	cout << "SampleRate: " << sampleRate << ", bitDepth: " << bitDepth << ", numSamples: " << numSamples << ", lenghtS: " << lengthInSeconds << '\n';		// Checking values
	cout << "Stark over " << '\n';																	// End of Checking Stark values
	cout << '\n';

	amplitudes.resize(new_size);																										// Resizing real data to nearest number power of 2
	std::vector<double> imag{ 0.0 };																									// Imag points
	imag.resize(new_size);																												// Resizing imag data to nearest number power of 2
	vector<complex<double>> cvec(amplitudes.size());																					// Input for fft
	vector<complex<double>> fft(amplitudes.size());																						// FFT output
	vector<complex<double>> ifft(amplitudes.size());																					// iFFT output

	transform(amplitudes.begin(), amplitudes.end(), imag.begin(), cvec.begin(), [](double da, double db)
	{
		return complex<double>(da, db);
	});

	//print("Original:", cvec.data(), log2N);
	cout << '\n';

	// Forward transform
	FFT(cvec.data(), fft.data(), log2N);
	//print("\nTransform:", fft.data(), log2N);

	// Generate an array, which contains the frequencies where the gaussian function will be generated. It contains 10 elements
	int gauss_freq[10] = { 31, 62, 125, 250, 500, 1000, 2000, 4000, 8000, 16000 };

	// Calculate the gaussian function
	int g = 0;
	double gaussian_function = 0;
	while (g < 10)
	{
		//gaussian_function += amplitude[i] * np.exp(-np.square(freq - gauss_freq[i]) / (2 * (gauss_freq[i] / 3)**2));
		g = g + 1;
	}

	//apply gaussian filter
	//for(int i = 0; i < N; i++ ) { ft[i] += ft[i] * gaussian_function; }

	iFFT(fft.data(), ifft.data(), log2N);
	//print("\nInvert:", ifft.data(), log2N);

	//-------------------------------------------------------------------------------------------------------
	//vector<complex<double>> *new_ifft;
	//std::vector<complex<double>> new_ifft = ifft.data;
	/*
	short int *value = new short int[samples_count];
	//Reading data
	for (int i = 0; i < samples_count; i++)
	{
		fread(&value[i], sample_size, 1, fin);
	*/
	std::vector<complex<double>> *new_ifft = ifft.data[sampleRate];

	for (int i = 0; i < wav.getNumSamplesPerChannel(); i++)
	{
		wav.samples[i] = &new_ifft;
	};
	//--------------------------------------------------------------------------------------------------------

	std::string outputFilePath = "quieter-audio-filer.wav"; // change this to somewhere useful for you
	wav.save(outputFilePath, AudioFileFormat::Wave);


	//WavOpen("input.wav", "output2.wav", 1, amplitudes, ifft);		// 1 for write, unified function
	cout << '\n';
	system("pause");

	return 0;
}

// Print values
void print(const char * prompt, Complex A[], int log2N)
{
	int N = 1 << log2N;
	cout << prompt << '\n' << fixed;
	for (int i = 0; i < N; i++) cout << A[i] << '\n';
}

// Fast Fourier Transform
void FFT(Complex f[], Complex ftilde[], int log2N)
{
	int N = 1 << log2N;

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
	}

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

// Finding nearest number power of 2
int resize(int samples)
{
	int nth_degree[21] = { 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824 };

	for (int i = 0; i < 21; i++)
	{
		if (samples < nth_degree[i])
		{
			return nth_degree[i];
			break;
		}
	}
}

// Wav file open and read
/*
-openwav(open_file, write_file, mode(1 for read, 2 for write - make it with if), vector)
*/
void WavOpen(const char* fileName, const char* writeName, int mode, vector<double> &V, vector<complex<double>> &iFFT)
{
	int read_or_write = mode;					// 1 - read wav data, 2 - write data to new wav 

	FILE *fin = fopen(fileName, "rb");

	//Read WAV header
	wav_header_t header;
	fread(&header, sizeof(header), 1, fin);

	//Print WAV header
	printf("WAV File Header read:\n");
	printf("File Type: %s\n", header.chunkID);
	printf("File Size: %ld\n", header.chunkSize);
	printf("WAV Marker: %s\n", header.format);
	printf("Format Name: %s\n", header.subchunk1ID);
	printf("Format Length: %ld\n", header.subchunk1Size);
	printf("Format Type: %hd\n", header.audioFormat);
	printf("Number of Channels: %hd\n", header.numChannels);
	printf("Sample Rate: %ld\n", header.sampleRate);
	printf("Sample Rate * Bits/Sample * Channels / 8: %ld\n", header.byteRate);
	printf("Bits per Sample * Channels / 8.1: %hd\n", header.blockAlign);
	printf("Bits per Sample: %hd\n", header.bitsPerSample);

	//skip wExtraFormatBytes & extra format bytes
	//fseek(f, header.chunkSize - 16, SEEK_CUR);

	//Reading file
	chunk_t chunk;
	printf("id\t" "size\n");

	//go to data chunk
	while (true)
	{
		fread(&chunk, sizeof(chunk), 1, fin);
		//printf("%c\t" "%li\n", chunk.ID[2], chunk.size);
		printf("%c%c%c%c\t" "%li\n", chunk.ID[0], chunk.ID[1], chunk.ID[2], chunk.ID[3], chunk.size);
		if (*(unsigned int *)&chunk.ID == 0x61746164)
			break;
		//skip chunk data bytes

		fseek(fin, chunk.size, SEEK_CUR);
	}

	//Number of samples
	int sample_size = header.bitsPerSample / 8;
	int samples_count = chunk.size * 8 / header.bitsPerSample;
	printf("Samples count = %i\n", samples_count);

	short int *value = new short int[samples_count];
	memset(value, 0, sizeof(short int) * samples_count);

	//Reading data
	for (int i = 0; i < samples_count; i++)
	{
		fread(&value[i], sample_size, 1, fin);
	}

	if (read_or_write == 0)
	{
		//Write data into the vector
		int amplitudes;
		for (int i = 0; i < samples_count; i++)
		{
			V.push_back(value[i]);
		}
		fclose(fin);
	}
	else if (read_or_write == 1)
	{
		printf("\nWriting new Wav started. \n");
		/*
		FILE* fp = fopen(writeName, "wb");

		if (!fp)
		{
		fprintf(stderr, "Couldn't open the file \"%s\"\n", writeName);
		//exit(0);
		system("pause");
		}
		*/

		ofstream NewWAV(writeName);									// Kiir wavba, nem lehet lejatszani
		for (int i = 0; i < samples_count; i++)						// Atirni a samples countot a new_size ertekre !!!!
		{
			value[i] = real(iFFT[i]);
			NewWAV << value[i];
		}

		printf("\nWriting new WAV ended successfully. \n");
		NewWAV.close();
		//fclose(fp);
	}
	else
	{
		printf("\nIncorrect use of function. \n		Use 0 for reading data from WAV \n		Use 1 for writing data to WAV \n");
	}
}