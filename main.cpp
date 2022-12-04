#include "AudioFile.h"
#include "fft.hpp"
#include "filter.h"

using namespace std;
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

int main()
{
    const std::string filePath = "guitar_mono.wav"; // get input file

    AudioFile<double> audioFile;              // audiFile (vector)
    vector<complex<double>> audioFileComplex; // audiofile converted to complex
                                              // so it can be used in in fft function
    vector<complex<double>> transform;
    size_t oldSize;
    size_t newSize;
    float sampleRate;

    bool loadedOK = audioFile.load(filePath); // read file

    if (!loadedOK)
        cout << "problem with reading the wav file";

    oldSize = audioFile.getNumSamplesPerChannel(); // save original size
    sampleRate = audioFile.getSampleRate();

    // resize to the closest 2 exponential, fill with 0's
    newSize = resize(audioFile.getNumSamplesPerChannel());
    audioFile.samples[0].resize(newSize);
    for (int i = oldSize; i < newSize; i++)
        audioFile.samples[0][i] = 0;

    // convert double to complex double
    // fill imaginary numbers with 0
    for (int i = 0; i < audioFile.getNumSamplesPerChannel(); i++)
        audioFileComplex.push_back(audioFile.samples[0][i]);

    transform = fft(audioFileComplex); // do the fft

    vector<double> freq(audioFile.getNumSamplesPerChannel());
    vector<double> gauss;
    vector<complex<double>> filteredTransform(transform);
    vector<complex<double>> ifftOutput;
    vector<double> filteredData(oldSize);

    // get  array with the list of frequencies
    freq = getFreq(audioFile.getNumSamplesPerChannel(), audioFile.getSampleRate());

    // generate the gaussian function
    gauss = generate_gaussian(freq);

    // apply the gaussian function to the signal
    for (int i = 0; i < freq.size(); i++)
        filteredTransform[i] += filteredTransform[i] * gauss[i];

    // normalize the signal so we dont have to much gain
    filteredTransform = normalizePower(filteredTransform, transform);

    AudioFile<double> outputAudio;

    // do the ifft
    ifftOutput = ifft(filteredTransform);

    // revert to the original size
    ifftOutput.resize(oldSize);

    // set to stereo
    outputAudio.setNumChannels(2);

    outputAudio.setNumSamplesPerChannel(ifftOutput.size());

    // set the sampleRate to the sample rate we read in the beginning
    outputAudio.setSampleRate(sampleRate);

    // copy the contents of the ifftouput array to the outputAudio object which can than be convert to a wav file
    for (int i = 0; i < ifftOutput.size(); i++)
    {
        for (int channel = 0; channel < outputAudio.getNumChannels(); channel++)
        {
            outputAudio.samples[channel][i] = ifftOutput[i].real();
        }
    }

    outputAudio.save("output.wav", AudioFileFormat::Wave);
}
