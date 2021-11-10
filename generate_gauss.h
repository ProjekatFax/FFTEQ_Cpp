#include <cmath>
#include <iostream>

/*param: 
		gauss - array with all values set to zero;
		gauss_freq - array with frequences where we generate our gaussian function
		gauss_amp - amplitude of a gauss function on a given freq
*/
void generate_gaussian(float* gauss, int* gauss_freq,float* gauss_amp){ 
	
	for(int i = 0; i<10; i++){
		for(int j = 0; j<100; j++)
			gauss[j] += gauss_amp[i] * exp(-(pow(j-gauss_freq[i],2) / (2*pow(gauss_freq[i]/3, 2))));
		
	}
}

