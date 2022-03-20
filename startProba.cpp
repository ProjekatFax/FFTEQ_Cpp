#include <complex>
#include <cstdio>
#include <cmath>
#include "string.h"
#include "memory.h"
#include "AudioFile.h"
#include "fft.hpp"
using namespace std;
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
    return samples;
}

int fftFreqSample(int sample){
    return sample/2 + 1;
}

int main(){
    const std::string filePath = "input.wav"; //get input file

    
    AudioFile<double> audioFile;  // audiFile (vector)
    vector<complex<double>> audioFileComplex; //audiofile converted to complex 
                                              //so it can be used in in fft function


    bool loadedOK = audioFile.load (filePath); //read file
    if(!loadedOK)
        cout<<"problem with reading the wav file";
    
    
    size_t newSize = resize(audioFile.getNumSamplesPerChannel());  
    audioFile.samples[0].resize(newSize);                   //resize to the closes 2 exponential, fill with 0's 
    int log2N = log2(newSize);	


    for(int i = 0; i<audioFile.getNumSamplesPerChannel(); i++)
        audioFileComplex.push_back(audioFile.samples[0][i]);   //convert double to complex double
                                                               //fill imaginary numbers with 0
    
    
    //int numberOfFrequencySamples = fftFreqSample(audioFile.getNumSamplesPerChannel()); 
    vector<complex<double>> audioFileFFT(audioFile.getNumSamplesPerChannel());
    
    fft(audioFileComplex.data(), audioFile.getNumSamplesPerChannel());
    
    
    
}





    

