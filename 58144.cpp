#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double euler(double y[], double x0, double (*function)(double, double, double), 
           int n, double h, double k, std::fstream &file);
double radioativeDecay(double k, double x, double y);
double rungeKutta(double y[], double x0, double (*function)(double, double, double),
                  int n, double h, double k, std::fstream &file);

int main(void){

int n=50;

double x0=0;
double h=0.5;
double y[n]={1,0};
double k=-2.3;
double erroEuler[n]={0};
double erroKutta[n]={0};

fstream decayEuler1;
fstream decayEuler2;
fstream decayEuler3;
decayEuler1.open("decayEuler1.dat", ios::out);
decayEuler2.open("decayEuler2.dat", ios::out);
decayEuler3.open("decayEuler3.dat", ios::out);

cout << "Exercicio 3.a)" << endl;
euler(y, x0, &radioativeDecay, n, h, k, decayEuler1);

h=0.7;
euler(y, x0, &radioativeDecay, n, h, k, decayEuler2);

h=1;
euler(y, x0, &radioativeDecay, n, h, k, decayEuler3);

cout << "Exercicio 3.b)" << endl;

fstream fErroEuler;
fstream fErroKutta;

fErroEuler.open("ErroEuler.dat", ios::out);
fErroKutta.open("ErroKutta.dat", ios::out);

int i = 0;
n=5;
for(double h=1; h>0.0078125; h = h/2){
    erroEuler[i] = euler(y, x0, &radioativeDecay, n, h, k, decayEuler3);
    erroKutta[i] = rungeKutta(y, x0, &radioativeDecay, n, h, k, decayEuler3);
    i++; 
}

for(i=0, h=1; i<7; i++){
    fErroEuler << h << "\t" << abs(erroEuler[i])-exp(-2.3*5) << endl;
    fErroKutta << h << "\t" << abs(erroKutta[i])-exp(-2.3*5) << endl;
    h/=2;
}


return 0;
}


double euler(double y[], double x0, double (*function)(double, double, double), 
int n, double h, double k, std::fstream &file){

for(int i=0; i<n; i++){
    y[i+1]=y[i]+h*(*function)(k, x0, y[i]);
    x0 += h;
    if(k!=5)
    file << y[i] << endl;
}

return y[n];
}

double radioativeDecay(double k, double x, double y){
    return k*y;
}

double rungeKutta(double y[], double x0, double (*function)(double, double, double),
int n, double h, double k, std::fstream &file){

for(int i=0; i<n; i++){
    double k1 = (*function)(k, x0, y[i]);
    double k2 = (*function)(k, x0+h/2, y[i]+k1*h/2);
    double k3 = (*function)(k, x0+h/2, y[i]+k2*h/2);
    double k4 = (*function)(k, x0+h, y[i]+k3*h);

    y[i+1]= y[i]+(h/6)*(k1+2*k2+2*k3+k4);
    x0=x0+h;

    if(k!=5)
    file << y[i] << endl;
}
return y[n];
}