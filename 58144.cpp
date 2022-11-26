#include <iostream>
#include <fstream>
using namespace std;

void euler(double y[], double x0, double (*function)(double, double, double), 
           int n, double h, double k, std::fstream &file);
double radioativeDecay(double k, double x, double y);


int main(void){

double x0=0;
double h=0.5;
int n=50;
double y[n]={1,0};

fstream Decay;
Decay.open("Decay.dat", ios::out);

euler(y, x0, &radioativeDecay, n, h, -2.3, Decay);




return 0;
}


void euler(double y[], double x0, double (*function)(double, double, double), 
           int n, double h, double k, std::fstream &file){

for(int i=0; i<n; i++){
    y[i+1]=y[i]+h*(*function)(k, x0, y[i]);
    x0 += h;
    file << y[i+1] << endl;
}
}

double radioativeDecay(double k, double x, double y){
    return k*y;
}
