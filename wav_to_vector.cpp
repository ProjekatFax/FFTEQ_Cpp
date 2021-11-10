//	https://github.com/tkaczenko/WavReader/blob/master/WavReader/WavReader.cpp
//  https://cplusplus.com/forum/general/265589/

#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <cstdint>
#include <stdio.h> 
#include <vector>
#include <iomanip>
#include <complex>
#include <cstdio>
#include <cmath>
#include "string.h"
#include "memory.h"

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
void WavReader(const char* fileName, vector<int> &V);
void print(const char * prompt, Complex A[], int log2N);
void FFT(Complex f[], Complex ftilde[], int log2N);
void iFFT(Complex ftilde[], Complex f[], int log2N);
int dataPoints();
int getFileSize(FILE* inFile);
//void add_to_array(int *&array, int size, int value)

int main(int argc, char* argv[])
{
	vector<int> amplitudes;

	WavReader("C:/Users/M�sz�ros Arnold/Documents/visual studio 2015/Projects/Project4/input.wav", amplitudes);			// Vraca data vrednosti u vektoru
	
	// For vector test
	int cnt = 0;
	for (int i = 0; i < amplitudes.size(); i++)
	{
		if (amplitudes[i] <= 32) 
		{
			cnt++;
		}
	}
	size_t size = amplitudes.size();
	cout << "Number of Vector elements: " << size << '\n';

	// Zero stuffing za sekvence koji imaju razlicit broj clanova od 2^n :    1, 2  ->  1, 2, 0, 0
	//add_to_array(data, D, 0)  - stavlja jednu nulu na mesto D+1. Resiti for petljom da dodaje onoliko 0 koliko treba

	//cout << "Number of rows in data.txt: " << countLines << " -- ovo je za debug funkcije fft";
	cout << '\n';

	//int fileLines = 1500000;				//---OVO JE const int N -> iz ovoga sledi log2N -- mozda

	//------------------------------------  Llink wav file to FFT
	const int N = 1498752, log2N = 4;			// Hardwired for testing, Values for testing the function - zameniti na dataPoints()

	Complex f[N] = amplitudes;
	Complex ft[N];		//array for FFTdata
	Complex n[N];		//array for IFFTdata

	print("Original:", f, log2N);
	cout << '\n';

	// Forward transform
	FFT(f, ft, log2N);   print("\nTransform:", ft, log2N);

	//freq = np.fft.rfftfreq(len(data), d = 1. / rate)


	//------------------------------------  Filter FFT_data

	// Generate an array, which contains the frequencies where the gaussian function will be generated. It contains 10 elements
	int gauss_freq[10] = { 31, 62, 125, 250, 500, 1000, 2000, 4000, 8000, 16000 };

	// Set preset defined in argument
	//amplitude = presets.set_presets(preset_arg)

	// Calculate the gaussian function
	int g = 0;
	double gaussian_function = 0;
	while (g < 10)
	{
		//gaussian_function += amplitude[i] * np.exp(-np.square(freq - gauss_freq[i]) / (2 * (gauss_freq[i] / 3)**2));
		g = g + 1;
	}

	//apply gaussian filter
	//FFT_data += FFT_data  * gaussian_function;

	//------------------------------------  Check FFT[].cpp == FFT_data.txt
	// Ovo samo debug

	//------------------------------------  ifft
	//iFFT(ft, f, log2N);   print("\nInvert:", f, log2N);
	iFFT(ft, n, log2N);   print("\nInvert:", n, log2N);
	
	//------------------------------------  Check ifft[].cpp == newdata.txt
	// Ovo isto ne treba ako su spojeni wav vrednosti sa IFFT
	
	// Isipis IFFT u wav


	cout << '\n';
	//fclose(wavFile);
	system("pause");

	return 0;
}

// Find the file size
int getFileSize(FILE* inFile)
{
	int fileSize = 0;
	fseek(inFile, 0, SEEK_END);

	fileSize = ftell(inFile);

	fseek(inFile, 0, SEEK_SET);
	return fileSize;
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

//	Count the lines in a file ( counting time points)
int dataPoints()
{
	FILE* fp;
	int countLines = 1;
	int i;
	int columns = 0;
	fp = fopen("data.txt", "r");
	if (fp == NULL)
	{
		cout << "Cant open file";
		fclose(fp);
		return 1;
	}
	else
	{
		while ((i = fgetc(fp)) != EOF)
		{
			if (i == '\n')
			{
				countLines++;
			}
		}

	}
	return countLines;
}

// Zero stuffing  - treba da prima vrednost N, da racuna razliku 2^n - N i to da doda na kraj novog reda
void add_to_array(int *&array, int size, int value) // add & to make "array" a reference
{
	int *newArr = new int[size + 1];
	memcpy(newArr, array, size * sizeof(int));
	delete[] array;
	array = newArr;
	array[size] = value;
} // Popraviti ovu funkciju

// Wav file open and read
void WavReader(const char* fileName, vector<int> &V)
{
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

	//Write data into the vector
	int amplitudes;
	for (int i = 0; i < samples_count; i++)
	{
		V.push_back(value[i]);
	}
	
}