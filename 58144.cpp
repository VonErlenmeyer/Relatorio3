#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

double metodoEuler(fstream &file, double x[], double y[], double h, double t0, int exercicio, bool File);
double metodoKutta(fstream &file, double x[], double y[], double h, double t0, int exercicio);


int main(void){

double x[500]={0};
double y[500]={0};
double h = 0.5;
double tInicial = 0;

fstream decayEuler1;
fstream decayEuler2;
fstream decayEuler3;
fstream trajetoriaEuler;
fstream trajetoriaKutta;
fstream erroEuler;
fstream erroKutta;

cout << "Exercicio 3. a)" << endl;

decayEuler1.open("DecayEuler1.dat", ios::out);
decayEuler2.open("DecayEuler2.dat", ios::out);
decayEuler3.open("DecayEuler3.dat", ios::out);

y[0] = 1;
metodoEuler(decayEuler1, x, y, 0.5, tInicial, 3, true);
metodoEuler(decayEuler2, x, y, 0.7, tInicial, 3, true);
metodoEuler(decayEuler3, x, y, 1, tInicial, 3, true);

cout << "Exercicio 3. b)" << endl;

erroEuler.open("ErroEuler.dat", ios::out);
erroKutta.open("ErroKutta.dat", ios::out);

for(double h=1; h>0.0078125; h/=2){
erroEuler << h << "\t" << fabs(metodoEuler(decayEuler3, x, y, h, tInicial, 3, false)-exp(-2.3*5)) << endl;
erroKutta << h << "\t" << fabs(metodoKutta(decayEuler3, x, y, h, tInicial, 3)-exp(-2.3*5)) << endl;
}

cout << "Exercicio 4" << endl;

y[0] = 0;
trajetoriaEuler.open("EulerTrajetoria.dat", ios::out);
trajetoriaKutta.open("KuttaTrajetoria.dat", ios::out);

metodoEuler(trajetoriaEuler, x, y, h, tInicial, 4, true);
metodoKutta(trajetoriaKutta, x, y, h, tInicial, 4);


}


double metodoEuler(fstream &file, double x[], double y[], double h, double t0, int exercicio, bool File){

int i=0;

if(exercicio == 3){
    for(i=0; i<=(5/h); i++){
    y[i+1]=y[i]+(-2.3*y[i]*h);
    if(File == true){
        file << t0 << "\t" << y[i] << endl;
        }
    t0 += h;
    }
}
else if (exercicio == 4){
    for(i=0; i<10; i++){
    y[i+1]=y[i]+(-9.81*t0+sqrt(200))*h;
    x[i+1]=x[i]+sqrt(200)*h;
    file << x[i] << "\t" << y[i] << endl;
    t0 += h;
    }
}

return y[i];
}

double metodoKutta(fstream &file, double x[], double y[], double h, double t0, int exercicio){

int i=0;
double yk1, yk2, yk3, yk4;
double xk1, xk2, xk3, xk4;

if(exercicio == 3){
    for(i=0; i<(5/h); i++){
    yk1 = (-2.3*(y[i]));
    yk2 = (-2.3*(y[i]+yk1*h/2));
    yk3 = (-2.3*(y[i]+yk2*h/2));
    yk4 = (-2.3*(y[i]+yk3*h));

    y[i+1] = y[i] + (h/6)*(yk1+2*yk2+2*yk3+yk4);

    }
}

else if (exercicio == 4){
    for(i=0;i<10; i++){
    yk1 = (-9.81*t0+sqrt(200));
    yk2 = (-9.81*(t0+h/2)+sqrt(200));
    yk3 = (-9.81*(t0+h/2)+sqrt(200));
    yk4 = (-9.81*(t0+h)+sqrt(200));

    xk1 = (sqrt(200));
    xk2 = (sqrt(200));
    xk3 = (sqrt(200));
    xk4 = (sqrt(200));

    x[i+1] = x[i] + (h/6)*(xk1+2*xk2+2*xk3+xk4);
    y[i+1] = y[i] + (h/6)*(yk1+2*yk2+2*yk3+yk4); 

    file << x[i] << "\t" << y[i] << endl;

    t0 += h;
    }
}

return y[i];
}


