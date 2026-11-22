#include <iostream> //general
#include <cstdlib>//for random numbers
#include <ctime> //for random numbers
#include <cmath> // math functions like exp, pow, etc
#include <random>// distributions etc
#include <algorithm> // for std::max, std::min, etc
#include <vector> // to use vectors, we need to include this library
#include <fstream> // for files
#include <iomanip> //for number of sf

//separate into order one odes and treat it like a vector. TAKES A BIT OF TIME.

std::vector<double> ODE(double x, std::vector<double> y){ //main function fn(x , vector)
    //now define the  vector, 
    std::vector<double> dy(4);
    dy[0]=y[1];
    dy[1]=y[2];
    dy[2]=y[3];
    dy[3]= (-50*y[2]-x*(y[0]*y[0]*y[0]))/10;
    return dy;
    
}
//we have ode now implement rk4. couldve done the past thing without vector but itll get too messsy

double rk4 (double x0, double xfina, double h){
    int N= std::round((xfina-x0)/h);
    std::vector<double> Yini={1,0,0,0};
    double x=x0;
    //we have to define a new vector that stores the intermediate steps. 
    for (int i=0; i<N; i++){
        //we now do the rk but must be careful, it has to be updated for each vector so we hqve to define it as vecotr
        //we here have to create the new auxiliary vector, we can not use Yini because that is for initial values and updated values ONLY. 
        //the auxiliary vector will store intermediate steps. very necessary. 
        //vectors need special care. of course there is the no vectors but define everything all time version, but this is cooler.
        std:: vector<double> intermediate(4);
        std::vector<double> k1 = ODE(x,Yini);
        //NOW ITS easy
        for (int j=0; j<4;j++) {
            intermediate[j]=Yini[j]+0.5*h*k1[j];
        }
        std:: vector<double> k2 =ODE (x+0.5*h, intermediate);// we compute using intermnediate
        for (int j=0; j<4;j++) {
            intermediate[j]=Yini[j]+0.5*h*k2[j];}
        std:: vector<double> k3 =ODE (x+0.5*h, intermediate);
        for (int j=0; j<4;j++) {
            intermediate[j]=Yini[j]+h*k3[j];}
        std::vector<double> k4 = ODE(x+h, intermediate);

        for (int j=0; j<4; j++) {
            Yini[j]= Yini[j]+h*(k1[j] + 2.0*k2[j]+ 2.0*k3[j] + k4[j])/6.0;
        }





        x=x+h;
    }
    return Yini[0];

}

int main() {
    double h=0.1;
    //can just approximate a random h or do a for loop that divides h by two each time until, say, 10 times and then compare
    //lets just do this for now
    std::cout <<std::setprecision(6);
    //try to do loop andsee values for different hs
    for (int k=0; k<10; k++){
      std::cout <<"Value for y(30) with six significant figures is "<< rk4(0,30,h)<<"." <<"For the value of h to be" << h <<std::endl;
      h= h/2.0; ///SEMICOLON!!!!!!!
    }

    
    return 0; //REMEMBER ; ALWAYS!!!
}
