#include <iostream> //general
#include <cstdlib>//for random numbers
#include <ctime> //for random numbers
#include <cmath> // math functions like exp, pow, etc
#include <random>// distributions etc
#include <algorithm> // for std::max, std::min, etc
#include <vector> // to use vectors, we need to include this library
#include <fstream> // for files

double t_initial = 0;
double t_final = 2.25;

double ode(double t, double x){
    return std::pow((t-1),2)*(x+1);
}
double exactsolution(double t) {
    return std::exp(t-std::pow(t,2)+(1.0/3.0)*std::pow(t,3))-1; //powers go with std too!!!. IF WE PUT DOUBLE WE HAVE TO PUT 1.0 OR WE CAN JSUT PUT 1??? put da zero
}

double RungeKutta2 (double t_initial, double t_final, double h , double x0) {

    int N = std::round((t_final - t_initial) / h); //dont put exact number just let it be smooth division, no pressure
    //initial conditions
    double x=x0;
    std::vector<double> x_rk2;// vectors to then define the function
    x_rk2.push_back(x0);
    double t= t_initial;
    std::vector<double> t_rk2;
    t_rk2.push_back(t_initial);

    for (int i=0; i<N; i++) {
        double k1=ode (t,x);
        double k2=ode(t+h, x+h*k1);
        x= x+h*(k1+k2)/2;
        t=t+h;
        x_rk2.push_back(x);
        t_rk2.push_back(t);

     


        }
        std::ofstream file("rk2_vector.csv"); //saves in csv for the graph in python
         file<< "t,x\n";// this is to put it in ordee
         for (int i=0; i<t_rk2.size(); i++) { // we create a loop that goes over the whole vectors
            file << t_rk2[i] << "," << x_rk2[i] << "\n"; } // t, x and line jump
        file.close(); //close. now what

   return x;
}




double RungeKutta4 (double t_initial, double t_final, double h , double x0) {
    int N = std::round((t_final - t_initial) / h);
    double x= x0;
    std::vector<double> x_rk4;
    x_rk4.push_back(x0);
    double t= t_initial;
    std::vector<double> t_rk4;
    t_rk4.push_back(t_initial);

    for (int i=0; i<N; i++) {
        double k1=ode (t,x); //this and k4 are the full step
        double k2=ode(t+(h/2), x+(h/2)*k1); // this and 3 are half step
        double k3=ode(t+(h/2), x+(h/2)*k2);
        double k4=ode(t+h, x+h*k3);
        x= x+h*(k1+2*k2+2*k3+k4)/6; // evaluate the weights, divide by the sum. 1+2+2+1=6
        t=t+h;
        
        x_rk4.push_back(x);
        t_rk4.push_back(t);
        
   }
   
   std::ofstream file("rk4_vector.csv");
         file<< "t,x\n";
         for (int i=0; i<t_rk4.size(); i++) {
            file << t_rk4[i] << "," << x_rk4[i] << "\n"; }
            
        file.close();
    
   return x;
}


void plot1 () { //this function is to get the first plot and the values of the solutions at h=0.25. 
    double h=0.25;
    double x0=0;
    double t_initial=0;
    double t_final=2.25; //the parameters. we want these for main

    double rk2=RungeKutta2(t_initial, t_final, h, x0);
    double rk4=RungeKutta4(t_initial, t_final, h, x0);
    double exact=exactsolution(t_final);
    
    std::vector<double> error_rk2;
    error_rk2.push_back(std::abs(rk2-exact));
    std::vector<double> error_rk4;
    error_rk4.push_back(std::abs(rk4-exact));
    std::vector<double> h_val;
    h_val.push_back(h);
    std::cout<<"Now, this is the First part"<<std::endl;
    std::cout << "THe exact solution is "<< exact << ". We can compare it with the RK2 Method whose solution give us" <<rk2<< " and can also conpare to RK4 whose solution is "<<rk4<<"."<< std::endl;
    std::cout << "The error for RK2 is "<< std::abs(rk2-exact) << " and the error for RK4 is "<< std::abs(rk4-exact) << "These are the final solutions, the file stores for each value of time, the progress. This file will appear in your VSC console as a csv file. "<<std::endl;}

void plot2() {//if we want the information on the error comparison for different h .  
    double x0=0;
    std::vector<double> error_rk2;
    
    std::vector<double> error_rk4;
    
    std::vector<double> h_val;
    std::cout<<"Sorry for the disordes, this is the Second part. It is the comparison of errors in RK2 compared to RK4 for different step sizes. "<<std::endl;
    for (double h=0.25; h>=0.0001; h=h/2) {
        double rk2=RungeKutta2(t_initial, t_final, h, x0);
        double rk4=RungeKutta4(t_initial, t_final, h, x0);
        double exact=exactsolution(t_final);
        error_rk2.push_back(std::abs(rk2-exact));
        error_rk4.push_back(std::abs(rk4-exact));
        h_val.push_back(h);
        
        
        std::cout << "For h= "<<h<<", the error for RK2 is "<< std::abs(rk2-exact) << " and the error for RK4 is "<< std::abs(rk4-exact) << std::endl;  
    }
    std::ofstream file("errortable.csv");
        file<< "h,error_rk2,error_rk4\n";
    for (int i=0; i<error_rk2.size();i++) {
            file <<h_val[i]<<","<<error_rk2[i]<<","<<error_rk4[i]<<"\n";

        }
        file.close();
}


//do the calls in main. i am sure there is a more effective way to do this
int main () {
    plot2();   //This is the second part, the errors for different step sizes. It mus be first because otherwise the graphs for the plot get solaped. 
    plot1();   //rk2_vector.csv  rk4_vector.csv. First part graphs. 
    
   
    
    return 0;
}