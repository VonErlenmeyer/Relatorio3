#include <iostream>
#include <fstream>
using namespace std;

int main(void){

double h=0.9;
double y[50]={0};
double x0 = 0;
y[0]=1;
int i;
double k1, k2, k3, k4;
fstream Test;
Test.open("Zeus.dat", ios::out);

for(i=0; i<50; i++){
    k1 = 5*y[i];
    k2 = 5*y[i]+k1*h/2;
    k3 = 5*y[i]+k2*h/2;
    k4 = 5*y[i]+k3*h;

    y[i+1]=y[i]+(h/6)*(k1+1*k2+2*k3+k4);
    x0 += h;
    Test << x0 << "\t" << y[i] << endl;
}

return 0;
}