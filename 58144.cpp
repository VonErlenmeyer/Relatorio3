#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

float trapezio(float a, float b, int n, float(*f)(float));
float simpson(float a, float b, int n, float(*f)(float));
float f(float x);
double trapezioArco(double F[], int n);
double velocidade(double m, double I);
double metodoEuler(fstream &file, double x[], double y[], double h, double t0, int exercicio, bool File);
double metodoKutta(fstream &file, double x[], double y[], double h, double t0, int exercicio);


int main(void){

int n;
float desvios;
float desviot;
float it;
float IT;

fstream Desviot;
Desviot.open("Desviot.dat", ios::out);

fstream Desvios;
Desvios.open("Desvios.dat", ios::out);

float Ia=123150;

for(n=2; n<=40;n=n+2){
    it = trapezio(4, 2, n, f);
    cout << "Valor do Integral: " << it << endl;
    desviot = fabs(it - Ia);
    cout << n << "\t" << desviot << endl;
    Desviot << n << "\t" << desviot << endl;
}

cout << "\n";

for(n=2; n<=40;n=n+2){
    IT = simpson(4, 2, n, f);
    cout << "Valor do Integral: " << IT << endl;
    desvios = fabs(IT - Ia);
    cout << n << "\t" << desvios << endl;
    Desvios << n << "\t" << desvios << endl;
    }

cout << "\nExercicio 2" << endl;

int k = 11;
double M = 0.075;
double F[11] = {0,37,71,104,134,161,185,207,225,239,250};
    
double I = trapezioArco(F,k);
cout << "I = " << I << endl;
velocidade(M,trapezioArco(F,k));

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

cout << "\nExercicio 3. a)" << endl;

decayEuler1.open("DecayEuler1.dat", ios::out);
decayEuler2.open("DecayEuler2.dat", ios::out);
decayEuler3.open("DecayEuler3.dat", ios::out);

y[0] = 1;
cout << "Para h=0.5:\n";
metodoEuler(decayEuler1, x, y, 0.5, tInicial, 3, true);
cout << "\nPara h=0.7\n";
metodoEuler(decayEuler2, x, y, 0.7, tInicial, 3, true);
cout << "\nPara h=1\n";
metodoEuler(decayEuler3, x, y, 1, tInicial, 3, true);

cout << "\nExercicio 3. b)" << endl;

erroEuler.open("ErroEuler.dat", ios::out);
erroKutta.open("ErroKutta.dat", ios::out);

for(double h=1; h>0.0078125; h/=2){
cout << "h = " << h << "\terro = " << fabs(metodoEuler(decayEuler3, x, y, h, tInicial, 3, false)-exp(-2.3*5)) << endl;
cout << "h = " << h << "\terro = " << fabs(metodoEuler(decayEuler3, x, y, h, tInicial, 3, false)-exp(-2.3*5)) << endl;
erroEuler << h << "\t" << fabs(metodoEuler(decayEuler3, x, y, h, tInicial, 3, false)-exp(-2.3*5)) << endl;
erroKutta << h << "\t" << fabs(metodoKutta(decayEuler3, x, y, h, tInicial, 3)-exp(-2.3*5)) << endl;
cout << "\n";
}

cout << "Exercicio 4" << endl;

y[0] = 0;
trajetoriaEuler.open("EulerTrajetoria.dat", ios::out);
trajetoriaKutta.open("KuttaTrajetoria.dat", ios::out);

cout << "Trajetoria, Euler (X,Y): " << endl;
metodoEuler(trajetoriaEuler, x, y, h, tInicial, 4, true);
cout << "\nTrajetoria, Kutta (X,Y): " << endl;
metodoKutta(trajetoriaKutta, x, y, h, tInicial, 4);

return 0;
}




float f(float x){
    return  -9800*(M_PI/3)*(32 - (x*x*(6 - x)));
}

float trapezio(float a, float b, int n, float(*f)(float)){

float h = (b-a)/n;
float soma_h = f(a)+f(b);
int i;
float Int;
    

for(i=1;i<n;i++){
soma_h=soma_h+2*f(a+i*h);
Int = (h/2)*soma_h;
}

return Int = (h/2)*soma_h;
}

float simpson(float a, float b, int n, float(*f)(float)){

float h = (b-a)/n;
float S_h = f(a)+f(b);
int i;
float Int;

for(i=1;i<n;i=i+2){
S_h=S_h+4*f(a+i*h);
Int = (h/3)*S_h;        
}

for(i=2;i<n-1;i=i+2){    
S_h=S_h+2*f(a+i*h);
Int = (h/3)*S_h;     
}

Int=(h/3)*(S_h);
        

return Int;
}

double velocidade(double m, double I){
double v = sqrt((I*2)/m);
cout << "v = " << v << endl;
return v;
}

double trapezioArco(double F[], int n){
    
double h = 0.05;
double S_h = F[n-1]+F[0];
    
for(int i = 1; i < n-1; i++){
    S_h=S_h+2*F[i];
}
double I = (h/2)*S_h;
//printf("I = %f\n",I);
return I; 
}

double metodoEuler(fstream &file, double x[], double y[], double h, double t0, int exercicio, bool File){

int i=0;

if(exercicio == 3){
    for(i=0; i<=(5/h); i++){
    y[i+1]=y[i]+(-2.3*y[i]*h);
    if(File == true){
        cout << t0 <<"\t" << y[i] << endl;
        file << t0 << "\t" << y[i] << endl;
        }
    t0 += h;
    }
}
else if (exercicio == 4){
double vy[300]={0};
    for(i=0; i<10; i++){
    vy[i+1]=vy[i]+(-9.81)*h;
    y[i+1]=y[i]+(vy[i]+sqrt(200))*h;
    x[i+1]=x[i]+sqrt(200)*h;
    cout << x[i] << "\t" << y[i] << endl;
    file << x[i] << "\t" << y[i] << endl;
    t0 += h;
    }
}

return y[i];
}

double metodoKutta(fstream &file, double x[], double y[], double h, double t0, int exercicio){

int i=0;
double vyk1, vyk2, vyk3, vyk4;
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
double vy[300]={0};
    for(i=0;i<10; i++){
    
    vyk1 = -9.81;
    vyk2 = -9.81;
    vyk3 = -9.81;
    vyk4 = -9.81;

    vy[i+1]=vy[i]+(h/6)*(vyk1*2+vyk2*2+vyk3+vyk4);

    yk1 = (vy[i]+sqrt(200));
    yk2 = (vy[i]+(-9.81*h/2)+sqrt(200));
    yk3 = (vy[i]+(-9.81*h/2)+sqrt(200));
    yk4 = (vy[i]+(h*(-9.81))+sqrt(200));

    xk1 = (sqrt(200));
    xk2 = (sqrt(200));
    xk3 = (sqrt(200));
    xk4 = (sqrt(200));

    x[i+1] = x[i] + (h/6)*(xk1+2*xk2+2*xk3+xk4);
    y[i+1] = y[i] + (h/6)*(yk1+2*yk2+2*yk3+yk4); 

    cout << x[i] << "\t" << y[i] << endl;
    file << x[i] << "\t" << y[i] << endl;

    t0 += h;
    }
}

return y[i];
}


