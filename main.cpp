#define SC_INCLUDE_FX
#include "AudioFile.h"
#include "fft.hpp"
#include "filter.h"
#include <deque>
#include <systemc>

using namespace std;

typedef sc_dt::sc_fix_fast num_t;
typedef deque<num_t> array_t;
typedef vector<double> orig_array_t;

void copy2fix(array_t &dest, const orig_array_t &src, int W, int F)
{
    for (size_t i = 0; i != src.size(); ++i)
    {
        num_t d(W, F);
        d = src[i];
        if (d.overflow_flag())
            std::cout << "Overflow in conversion.\n";
        dest.push_back(d);
    }
}

bool passCheck(const orig_array_t &gold, const orig_array_t &sys,
               double delta)
{
    for (size_t i = 0; i != gold.size(); ++i)
    {
        if (std::abs(gold[i] - sys[i]) > delta)
            return false;
    }
    return true;
}

vector<double> readFromFile(string fileName){
    // Vector to hold the doubles read from the file
    std::vector<double> numbers;

    // Open the text file in read mode
    std::ifstream infile("output.txt");

    // Check if the file stream is open
    if (!infile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
    }

    std::string line;
    // Read the file line by line
    while (std::getline(infile, line)) {
        // Convert line to a double
        std::istringstream iss(line);
        double number;
        if (iss >> number) {
            // Add the number to the vector
            numbers.push_back(number);
        } else {
            std::cerr << "Failed to convert line to double: " << line << std::endl;
        }
    }

    // Close the file
    infile.close();

    // Output the numbers to verify they were read correctly
    std::cout << "Numbers read from file:" << std::endl;
    for (double num : numbers) {
        std::cout << num << std::endl;
    }

    return numbers;
}

int main(int argc, char **argv)
{
    array_t input;
    uint8_t preset;
    string arg1 = argv[1]; // get input file
    string arg2;
    if (arg1 == "-h")
    {
        cout << "./filename 'wavFileNam' 'presets' " << endl;
        cout << "presets can be: bass, mid, trebble" << endl;
        cout << "example: ./program someFile.wav bass" << endl;
        return 0;
    }

    if (argc > 2)
        arg2 = argv[2];

    if (arg2 == "bass")
        preset = BASS;
    else if (arg2 == "mid")
        preset = MID;
    else if (arg2 == "treble")
        preset = TREBLE;
    else
    {
        cout << "Argument for presets incorrect" << endl;
        cout << "call with -h argument: ./filename -h" << endl;
        return 1;
    }
    AudioFile<double> audioFile;              // audiFile (vector)
    vector<complex<double>> audioFileComplex; // audiofile converted to complex
                                              // so it can be used in in fft function
    vector<complex<double>> transform;
    size_t oldSize;
    size_t newSize;
    float sampleRate;

    bool loadedOK = audioFile.load(arg1); // read file

    if (!loadedOK)
    {
        cout << "problem with reading the wav file";
        cout << "call with -h argument: ./filename -h" << endl;

        return 1;
    }   


    oldSize = audioFile.getNumSamplesPerChannel(); // save original size
    sampleRate = audioFile.getSampleRate();

    // resize to the closest 2 exponential, fill with 0's
    newSize = resize(audioFile.getNumSamplesPerChannel());
    audioFile.samples[0].resize(newSize);
    for (int i = oldSize; i < newSize; i++)
        audioFile.samples[0][i] = 0;
    bool pass = false;
    int W = 5;
    const int F = 2;
    const double error_d = 1e-3;
    vector<double> gold = readFromFile(arg2 + ".txt");



    do
    {
        copy2fix(input, audioFile.samples[0], W, F);

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
        gauss = generate_gaussian(freq, preset);

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

        //outputAudio.save("output.wav", AudioFileFormat::Wave);
        
        pass = passCheck(gold, sys, error_d);
        W++;
        input.clear();
    } while (pass == false);

    return 0;
}
