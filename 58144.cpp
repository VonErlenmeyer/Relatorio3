#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double euler(double y[], double x0, double (*function)(double, double), 
             double h, std::fstream &file, bool File);
double radioativeDecay(double x, double y);
double rungeKutta(double y[], double x0, double (*function)(double, double),
                  double h, std::fstream &file, bool File);

int main(void){

double x0=0;
double h=0.5;
double y[100]={1,0};
bool File=true;


fstream decayEuler1;
fstream decayEuler2;
fstream decayEuler3;
decayEuler1.open("decayEuler1.dat", ios::out);
decayEuler2.open("decayEuler2.dat", ios::out);
decayEuler3.open("decayEuler3.dat", ios::out);

cout << "Exercicio 3.a)" << endl;
euler(y, x0, &radioativeDecay, h, decayEuler1, File);

h=0.7;
euler(y, x0, &radioativeDecay, h, decayEuler2, File);

h=1;
euler(y, x0, &radioativeDecay, h, decayEuler3, File);

cout << "Exercicio 3.b)" << endl;

fstream fErroEuler;
fstream fErroKutta;

fErroEuler.open("ErroEuler.dat", ios::out);
fErroKutta.open("ErroKutta.dat", ios::out);

File = false;
for(double h=1; h>0.0078125; h = h/2){
fErroEuler << h << "\t" << abs(abs(euler(y, x0, &radioativeDecay, h, decayEuler3, File))-exp(-2.3*5)) << endl;
fErroKutta << h << "\t" << abs(abs(rungeKutta(y, x0, &radioativeDecay, h, decayEuler3, File))-exp(-2.3*5)) << endl;
}




return 0;
}


double euler(double y[], double x0, double (*function)(double, double), 
  double h, std::fstream &file, bool File){

int i=0;
for(; i<=(5/h); i++){
    y[i+1]=y[i]+h*(*function)(x0, y[i]);
    if(File == true){
    file << x0 << "\t" << y[i] << endl;
    }
    x0 += h;
}

return y[i-1];
}

double radioativeDecay(double x, double y){
    return -2.3*y;
}

double rungeKutta(double y[], double x0, double (*function)(double, double),
  double h, std::fstream &file, bool File){

int i=0;
for(; i<(5/h); i++){
    double k1 = (*function)(x0, y[i]);
    double k2 = (*function)(x0+h/2, y[i]+k1*h/2);
    double k3 = (*function)(x0+h/2, y[i]+k2*h/2);
    double k4 = (*function)(x0+h, y[i]+k3*h);

    y[i+1]= y[i]+(h/6)*(k1+2*k2+2*k3+k4);
    x0=x0+h;

    if(File == true)
    file << y[i] << endl;
}
return y[i-1];
}