#include <iostream>
#include <algorithm> // For std::transform
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
#include <math.h>
#include "string.h"
#include "memory.h"
#include "AudioFile.h"

using namespace std;

void writeSineWaveToAudioFile()
{
      AudioFile<float> a;
      a.setNumChannels (1);
      a.setNumSamplesPerChannel(1000000);

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

int main(){
    writeSineWaveToAudioFile();

    return 0;
}