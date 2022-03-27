#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

/*param: 
		gauss - array with all values set to zero;
		gauss_freq - array with frequences where we generate our gaussian function
		gauss_amp - amplitude of a gauss function on a given freq
*/
vector<int> gauss_freq {31,62,125,250,500,1000,2000,4000,8000,16000};

vector<int> gauss_amp {6,4,2,0,0,0,0,0,0,0};

void generate_gaussian(double  gauss[],int n, double freq[]){ 
	
	for(int i : gauss_freq){
		for(int j = 0; j<n; j++)
			gauss[j] += gauss_amp[i] * exp(-(pow(freq[j]-gauss_freq[i],2) / (2*pow(gauss_freq[i]/3, 2))));
	}
}


