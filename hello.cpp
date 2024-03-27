// Code your design here.
// Uncomment the next line for SystemC modules.
#define SC_INCLUDE_FX


#include <systemc>
#include <math.h>
#include <iostream>


int factorial(int n){
  	int result = 1;
  	
    while(n>0){
		result = result*n;
      	n--;
    }
  
  	return result;
}


double calculate_sine_float(float x,int k){
	int n;
  	double result = 0;	
  
  	for( n = 0; n<=k; n++){
  		result += pow(-1,n)/factorial(2*n+1)*pow(x,(2*n+1));
    }
  	
  	return result;
}

sc_dt::sc_fixed<70,58,sc_dt::SC_RND_CONV,sc_dt::SC_SAT,5>  calculate_sine_fixed(float x,int k){
	int n;
  	sc_dt::sc_fixed<60,50,sc_dt::SC_RND_CONV,sc_dt::SC_SAT,5> result = 0;	
  
  	for( n = 0; n<=k; n++){
  		result += pow(-1,n)/factorial(2*n+1)*pow(x,(2*n+1));
      	std::cout<<"loop number:"<< n << " " << result << std::endl;
    }
  	
  	return result;
}



int sc_main(int argc, char* argv[]){
	int n = 8;
  	float x_float = 0.5;
  	
  	double sine = calculate_sine_float(x_float,n);
    sc_dt::sc_fixed<70,58,sc_dt::SC_RND_CONV,sc_dt::SC_SAT,5> sine_fixed = calculate_sine_fixed(x_float,n);
  	std::cout <<"sine float: "<< sine <<std::endl;
  	std::cout <<"sine fixed: "<< sine_fixed <<std::endl;

	return 0;
}