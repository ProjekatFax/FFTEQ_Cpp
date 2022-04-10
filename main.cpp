#include <complex>
#include <cstdio>
#include <cmath>
#include "string.h"
#include "memory.h"
#include "AudioFile.h"
#include "fft.hpp"
#include "generate_gauss.h"

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


int main(){
    const std::string filePath = "sine-wave.wav"; //get input file

    
    AudioFile<double> audioFile;  // audiFile (vector)
    vector<complex<double>> audioFileComplex; //audiofile converted to complex 
                                              //so it can be used in in fft function

    AudioFile<double> outputAudio;

    

    bool loadedOK = audioFile.load (filePath); //read file
    if(!loadedOK)
        cout<<"problem with reading the wav file";
    
    size_t oldSize = audioFile.getNumSamplesPerChannel();
    
    size_t newSize = resize(audioFile.getNumSamplesPerChannel());  
    audioFile.samples[0].resize(newSize);                   //resize to the closes 2 exponential, fill with 0's 
    for(int i = oldSize; i<newSize; i++)
        audioFile.samples[0][i];

    int log2N = log2(newSize);	

    for(int i = 0; i<audioFile.getNumSamplesPerChannel(); i++)
        audioFileComplex.push_back(audioFile.samples[0][i]);   //convert double to complex double
                                                               //fill imaginary numbers with 0
    
    //int numberOfFrequencySamples = fftFreqSample(audioFile.getNumSamplesPerChannel()); 
    vector<complex<double>> audioFileFFT(audioFile.getNumSamplesPerChannel());
    
    fft(audioFileComplex.data(), audioFile.getNumSamplesPerChannel());

    
    vector<double> freq(audioFile.getNumSamplesPerChannel());
    vector<double> gauss(audioFile.getNumSamplesPerChannel());
    vector<complex<double>> filtered(audioFile.getNumSamplesPerChannel());
    vector<double> filteredData(audioFile.getNumSamplesPerChannel());
    vector<complex<double>> ifftOutput(audioFile.getNumSamplesPerChannel());



    getFreq(freq.data(), audioFile.getNumSamplesPerChannel(),audioFile.getSampleRate());

    generate_gaussian(gauss.data(), audioFile.getNumSamplesPerChannel(),freq.data());

    filtered = audioFileComplex;
    cout << "filtered at 0 before gauss : " <<filtered[0] << endl;

    
    //for(int i = 1; i<audioFile.getNumSamplesPerChannel(); i++)
    //    filtered[i] += filtered[i]*gauss[i];

    filtered[0] = 0;
    
    ifft(filtered.data(),ifftOutput.data(), audioFile.getNumSamplesPerChannel());

    cout << "ifft size" << filtered.size()<<endl;

    filtered.resize(oldSize);

    for(int i = 0; i<=oldSize;i++){
        filteredData[i] = filtered[i].real();
    }

    outputAudio.setNumChannels (1);
    outputAudio.setNumSamplesPerChannel (oldSize);


    for (int i = 0; i < outputAudio.getNumSamplesPerChannel(); i++)
        {
            for (int channel = 0; channel < outputAudio.getNumChannels(); channel++)
            {
                outputAudio.samples[channel][i] = filteredData[i];
            }
        }

    outputAudio.save("output.wav", AudioFileFormat::Wave);
}

void writeSineWaveToAudioFile()
{
      AudioFile<float> a;
      a.setNumChannels (1);
      a.setNumSamplesPerChannel(10000);

      const float sampleRate = 44800.f;
      const float freq = 440.f;

      for (int i = 0; i < a.getNumSamplesPerChannel(); i++)
      {
          for (int channel0 = 0; channel0 < a.getNumChannels(); channel0++){
            a.samples[channel0][i] = sin((static_cast<float> (i) / sampleRate) * freq * 2.f * M_PI);
          }
      }

       string filePath = "sine-wave.wav"; // change this to somewhere useful for you
        a.save ("sine-wave.wav", AudioFileFormat::Wave);

}

