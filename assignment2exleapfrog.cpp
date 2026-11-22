//hi
#include <iostream> //general
#include <cstdlib>//for random numbers
#include <ctime> //for random numbers
#include <cmath> // math functions like exp, pow, etc
#include <random>// distributions etc
#include <algorithm> // for std::max, std::min, etc
#include <vector> // to use vectors, we need to include this library
#include <fstream> // for files
#include <iomanip> //for precision 

struct planetinfo{
    std:: vector<double> x0;// xaxis
    std::vector<double> x1; //yaxis
    std:: vector<double> v0; // vx componen
    std:: vector<double> v1; // vy component DONT FORGET ;
    std:: vector<double> mass;
    std::vector<double> acceleration0; // axaxis
    std::vector<double> acceleration1;//yaxis

};// ; after struct ALWAYS, struct acts.
// planetinfo planet(){}would be a function but we re definign the variables for each planet in vector

void accelerations (planetinfo& plan) {// wont give anything just modifies variabes so its a void.

 for (int i=0; i<plan.mass.size(); i++){
    double acc0=0;
    double acc1=0;//each planet i 
    for (int j=0; j<plan.mass.size(); j++){///force fron other planets on each i
        if(i==j){//comparing i to an other planet j and get the force. 
            continue; //continue means skip that
        }
        //dont really need an else
        //rij is square root of change x squared plus change y squared
        double diffx= plan.x0[j]-plan.x0[i]; //for x
        double diffy= plan.x1[j]-plan.x1[i];//for y
        //acceleration is Gmj times r vector over r norm cube..
        //separare in two components. r vector for one id diffx for the other diffy. the r norm cube is same. 
        //G is one G=1
        double G =1;
        double rnormcube =std::pow(diffx*diffx+diffy*diffy, 1.5);//update accelerations
        acc0+= G*(plan.mass[j]*diffx)/rnormcube;
        acc1+=G*(plan.mass[j]*diffy)/rnormcube;

    }
    //get the updage
    plan.acceleration0[i]=acc0;
    plan.acceleration1[i]=acc1;

 }

}

void leapfrog (planetinfo& plan, double h){
    int numberplanet= plan.mass.size();
    for (int i = 0; i < numberplanet; i++) { //vels with the first acce
        plan.v0[i] += 0.5*h*plan.acceleration0[i];
        plan.v1[i] += 0.5*h*plan.acceleration1[i];
    }

    for (int i=0; i<numberplanet; i++){
        plan.x0[i] += h*plan.v0[i];//here we move
        plan.x1[i] += h*plan.v1[i];
    }
    //UPDATE ACCELERATIOONNN!!!
    accelerations(plan);//update acc because well, acc changes w positiron
    for (int i=0; i<numberplanet; i++){
        plan.v0[i] += 0.5*h*plan.acceleration0[i];//new half velocity now does loop
        plan.v1[i] += 0.5*h*plan.acceleration1[i];
    }

}




int main(){
planetinfo planet; //define like variable. be careful with this.
planet.x0= {-0.60, -0.50, 0.50, 0.50};
planet.x1={-0.20, 0.10, 0.10, 0.40};
planet.v0= {1.86, -0.84, -0.44, 1.15};
planet.v1= {0.70, 0.65,-1.50, -1.60};
planet.mass= {0.8, 2.2, 0.9, 0.4};
planet.acceleration0={0,0,0,0}; //compute now or re define later?
planet.acceleration1={0,0,0,0};
accelerations(planet);
for (int i = 0; i < 4; i++){//initiql acc
    std::cout << i << ": " << planet.acceleration0[i] << "  " << planet.acceleration1[i] << "\n";
}

double tfinal=4.5;
double Nsteps=100000;
double h= tfinal/Nsteps;
for (int k = 0; k < Nsteps; k++){
    
    leapfrog(planet,h);
    
}


for (int i=0;i<4;i++){
    std::cout <<"The planet "<<i<< "Has a position in the y-axis of "<<planet.x1[i]<<"and in the x-axis of  "<<planet.x0[i]<<"\n";
}

return 0;
}
