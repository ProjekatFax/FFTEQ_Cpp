#include <complex>
#include <iostream>
#include <vector>
using namespace std;

void magija(complex<int> *c,int size){
    c[0] = 2;
}


int main(){
    std::vector<complex<int>> c; 
    
    c.push_back(1);
    c.push_back(1i);
    complex<int>* cPointer = c.data();
    magija(cPointer, c.size());
    for(int i = 0; i<c.size();i++)
        cout<<"real: "<< c[i].real()<<" imag: " << c[i].imag() <<endl;
}

