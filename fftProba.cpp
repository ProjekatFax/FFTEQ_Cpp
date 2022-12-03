#include <complex>
#include <cstdio>
#include <cmath>
#include "string.h"
#include "memory.h"
#include "fft.hpp"

using namespace std;

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
    const int N = 8, log2N = 3;   
    vector<complex<double>> f{ 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 1.0, -1.0 };
    vector<complex<double>> ft;
    vector<complex<double>> output;

    cout<<"original: "<<endl;
    for(int i = 0; i < f.size() ; i++){
        cout<<f[i]<<endl;
    }
    cout<<endl;

    ft =  fft(f);

    cout<<"ft size: " <<ft.size()<<endl;

    for(int i = 0; i < ft.size() ; i++){
        cout<<ft[i]<<endl;
    }
    cout<<endl;

    output = ifft(ft);
    
    for(int i = 0; i < output.size() ; i++){
        output[i].imag(0);
        cout<<output[i]<<endl;
    }
}