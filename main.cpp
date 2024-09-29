#define SC_INCLUDE_FX
#include <systemc>
#include "AudioFile.h"
#include "fft.hpp"
#include "filter.h"
#include <deque>

using namespace std;

typedef sc_dt::sc_fix_fast num_t;
typedef deque<num_t> array_t;
typedef vector<double> orig_array_t;

void copy2fix(orig_array_t& dest, const orig_array_t& src, int width, int fraction)
{

    for (size_t i = 0; i != src.size(); ++i)
    {
        num_t delta(width, fraction);
        
        delta = src[i];
        if(i%1500000==0)
        cout << "vrednost Gold "<<src[i] << " vrednost crippled: " << delta << endl;
        if (delta.overflow_flag())
            std::cout << "Overflow in conversion.\n";
        dest.push_back(delta);
    }
}

bool passCheck(const orig_array_t gold, const orig_array_t sys,
               double delta)
{
    //cout<< "GOLD SIZE: " << gold.size() << endl;
    for (size_t i = 0; i != gold.size(); ++i)
    {

        double sum = std::abs(std::abs(gold[i]) - std::abs(sys[i]));
        if(i%1500000==0)
          cout<< "value of error: " << sum << endl;
        if ( sum > delta){
            cout << "error place: " << i << endl;
            cout << "value of gold[i]: " << std::abs(gold[i]) <<" - " << "sys[i]: " <<  std::abs(sys[i]) << endl;
            cout << "value of error: " << sum << endl;
            
            return false;
        }
    }
    cout << "passChech Success" << endl;
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
    return numbers;
}

void print_help() {

    cout << "./filename 'wavFileNam' 'presets' " << endl;
    cout << "presets can be: bass, mid, trebble" << endl;
    cout << "example: ./program someFile.wav bass" << endl;
}

void print_bad_preset_name(){
    cout << "Argument for presets incorrect" << endl;
    cout << "call with -h argument: ./filename -h" << endl;
}

bool load_audio_file(string filename, AudioFile<double> audioFile) {
    cout << "Loading audio file " << filename << " ..." << endl;
    return audioFile.load(filename); // read file
}

int sc_main(int argc, char **argv)
{
    AudioFile<double> audioFile;              // audiFile (vector)
    vector<complex<double>> audioFileComplex; // audiofile converted to complex
                                                // so it can be used in in fft function
    vector<complex<double>> transform;
    size_t oldSize;
    size_t newSize;
    float sampleRate;

    string filename = argv[1]; // file name
    if (filename == "-h")
    {
        print_help();
        return 0;
    }
    string preset_name;

    if (argc > 2)
        preset_name = argv[2];

    uint8_t preset;

    if (preset_name == "bass")
        preset = BASS;
    else if (preset_name == "mid")
        preset = MID;
    else if (preset_name == "treble")
        preset = TREBLE;
    else
    {
        print_bad_preset_name();
        return 1;
    }

    cout << "Loading audio file " << filename << " ..." << endl;
    bool isOk = audioFile.load(filename); // read file

    // if (!load_audio_file(filename, *audioFile))
    if (!isOk)
    {
        cout << "problem with  reading the wav file ";
        cout << "call with -h argument: ./filename -h" << endl;
        return 1;
    }  

    // resize_matrix // append zeros?

    oldSize = audioFile.getNumSamplesPerChannel(); // save original size
    sampleRate = audioFile.getSampleRate();

    //cout << "oldSize" << oldSize << endl;
    //cout << "sampleRate" << sampleRate << endl;

    // cout << "problem with  reading the wav file";
    // cout << "call with -h argument: ./filename -h" << endl;

    // resize to the closest 2 exponential, fill with 0's
    newSize = resize(audioFile.getNumSamplesPerChannel());
    cout << "newSize" << newSize << endl;
    audioFile.samples[0].resize(newSize);
    for (int i = oldSize; i < newSize; i++)
        audioFile.samples[0][i] = 0;

    // 
    
    bool pass = false;
    int width = 64;
    int fraction = 3;
    const double error_delta = 1e-3;

    orig_array_t sys;
    orig_array_t gold;

    int iteration_number = 0;

    orig_array_t input;
    do
    {
        cout << "Iteration" << iteration_number << endl;
        //cout << "Width" << width << endl;
        // ifftOutput.clear();
        input.clear();
        sys.clear();
        audioFileComplex.clear();

        copy2fix(input, audioFile.samples[0], width, fraction);  //we get the crippled version
        
        // convert double to complex double
        // fill imaginary numbers with 0
        //cout << "audioFile.getNumSamplesPerChannel() " << audioFile.getNumSamplesPerChannel() << endl;

        //cout << "Start audioFileComplex" << endl; 
        for (int i = 0; i < audioFile.getNumSamplesPerChannel(); i++){
            audioFileComplex.push_back(input[i]); 
        }
        // cout << "End audioFileComplex" << endl;
        
        
        // cout << "Start fft" << endl;
        transform = fft(audioFileComplex); // do the fft
        vector<double> freq(audioFile.getNumSamplesPerChannel());
        vector<double> gauss;
        vector<complex<double>> filteredTransform(transform);
        vector<complex<double>> ifftOutput;
        vector<double> filteredData(oldSize);
        //cout << "c1" << endl;

        // get  array with the list of frequencies
        freq = getFreq(audioFile.getNumSamplesPerChannel(), audioFile.getSampleRate());

        //cout << "c2" << endl; 

        // generate the gaussian function
        gauss = generate_gaussian(freq, preset);

        //cout << "c3" << endl;

        // apply the gaussian function to the signal
        for (int i = 0; i < freq.size(); i++)
            filteredTransform[i] += filteredTransform[i] * gauss[i];

        //cout << "c4" << endl;
        
        // normalize the signal so we dont have too much gain
        filteredTransform = normalizePower(filteredTransform, transform);

        //cout << "c5" << endl;

        AudioFile<double> outputAudio;

        // do the ifft
        ifftOutput = ifft(filteredTransform);

        //cout << "c6" << endl;

        // revert to the original size
        ifftOutput.resize(oldSize);

        //cout << "c7" << endl;

        // set to stereo
        outputAudio.setNumChannels(2);

        outputAudio.setNumSamplesPerChannel(ifftOutput.size());

        // set the sampleRate to the sample rate we read in the beginning
        outputAudio.setSampleRate(sampleRate);

        // copy the contents of the ifftouput array to the outputAudio object which can than be convert to a wav file
        //cout << "Start copy contents" << endl; 
        for (int i = 0; i < ifftOutput.size(); i++)
        {
            for (int channel = 0; channel < outputAudio.getNumChannels(); channel++)
            {    
                sys.push_back(ifftOutput[i].real());
                outputAudio.samples[channel][i] = ifftOutput[i].real();
            }
        } 
        cout << "End copy contents" << endl;

        

        if(iteration_number != 0){
            pass = passCheck(gold, sys, error_delta);
            // newSize = resize(audioFile.getNumSamplesPerChannel());
            // audioFile.samples[0].resize(newSize);
        }else{        
            cout << "-->GOLD<--" << endl; 
            // cout <<  passCheck(sys, sys, error_delta) <<  endl; 
            gold = sys;
            width =1;
        }
        width++;

        //For saving each iteration of the width into a .wav file
        //std::stringstream ss;
        //ss << "output_" << width << "_" << preset_name << ".wav";
        //std::string outputFilename = ss.str();

        /* if (width == 3) {
           outputAudio.save("output12.wav", AudioFileFormat::Wave);
        } else if (width == 5) {
            outputAudio.save("output13.wav", AudioFileFormat::Wave); 
        } else if (width == 10) {
            outputAudio.save("output14.wav", AudioFileFormat::Wave);
        } else if (width == 15) {
            outputAudio.save("output15.wav", AudioFileFormat::Wave);
        } */

        //outputAudio.save(outputFilename, AudioFileFormat::Wave);
        
        iteration_number++;

    } while (pass == false);

    // outputAudio.save("outputwwww.wav", AudioFileFormat::Wave)

    cout << "width: " << width<< endl;
    cout << "Fraction: " << fraction << endl;

    return 0;
}
