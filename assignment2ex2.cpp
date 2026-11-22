#include <iostream> //general
#include <cstdlib>//for random numbers
#include <ctime> //for random numbers
#include <cmath> // math functions like exp, pow, etc
#include <random>// distributions etc
#include <algorithm> // for std::max, std::min, etc
#include <vector> // to use vectors, we need to include this library
#include <fstream> // for files
#include <iomanip> //for precision lol

double t_initial = 0;
double t_final = 15;

double eqn1 (double v) { // v
    return v;
}

double eqn2 (double x, double t) {//x''
    return -(t*x*(t+3))/(2+(t*t)*(x*x));
}
struct allsolns {//define this type of thing, rhe solutions. struct works well to group properties.
    std::vector<double> tva;
    std::vector<double> xva;
    std::vector<double> vva;
};// has this here
allsolns rungekutta (double x0, double t_initial, double t_final, double h, double v0) {
    int N = std::round((t_final - t_initial) / h); //return and name of variable must be same. this struct thing is very interesting. keep using. 
                   
    double x = x0; 
    double v = v0; 
    double t = t_initial;
    allsolns sol; //like double but now is allsolns, like a type of variable, and sol the name
    

    sol.xva.push_back(x0);
    sol.vva.push_back(v0); 
    sol.tva.push_back(t_initial);
    

    for (int i = 0; i < N; i++) {
        double kx1 = eqn1(v);
        double kv1 = eqn2(x, t);

        double kx2 = eqn1(v + h/2*kv1);
        double kv2 = eqn2(x + h/2*kx1, t + h/2);

        double kx3 = eqn1(v + h/2*kv2);
        double kv3 = eqn2(x + h/2*kx2, t + h/2);

        double kx4 = eqn1(v + h*kv3);
        double kv4 = eqn2(x + h*kx3, t + h);

        x += h/6 * (kx1 + 2*kx2 + 2*kx3 + kx4);
        v += h/6 * (kv1 + 2*kv2 + 2*kv3 + kv4);
        t += h;

        sol.xva.push_back(x);
        sol.vva.push_back(v);
        sol.tva.push_back(t);

    }

    
    return sol;//can return sol allsols bc insted of defining func as double i deinfed as allsols, and defined all sols by struct. 

}



//lets do shooting method

double errorfunction (double v0, double h) {
    return rungekutta(1, 0, 15, h, v0).xva.back()-1;// want to find its roots
}// cn call a struct in double

double bisectionmethod (double a, double b, double h){
 double fa= errorfunction (a,h);
 double fb= errorfunction (b,h);
//check for change in sign. if there IS the prod will be negative otherwise orod will be posirive

    if (fa*fb>0) {
        std::cout<<"No change of sign between"<<a <<"and"<<b<<"\n!";
        return NAN;
    } //nan means nothin

    for (int i=0; i<100; i++) {
     double m = 0.5*(a+b);
     double fm=errorfunction(m,h); 
        if (fa*fm<=0){
        //negaative theere is root
        b=m;
        fb=fm;
    }
    else {
    //we looking for root here
    a=m;
    fa=fm;
    }
    }
return 0.5*(a+b);//not return m m is defined in the for, return w the updated values of a,b.
}

std::vector<double> allroots() { // roots and then we use them in rk to get the different solutios
    double vmin=-100;
    double dv= 0.005;
    double vfinal = 100;
    double h=0.01; //for rk
    double steps = (std::abs(vfinal)+std::abs(vmin))/dv;
    double vprev=vmin;
    double errfprev= errorfunction(vprev,h);
    std::vector<double> roots;
    for (int j=0; j<=steps; j++){
        double ve=vmin+j*dv;//mult by j to get a loop othwerwise the dv value is always the same
        double errfnow= errorfunction(ve,h);
        if (errfprev*errfnow<0) {
            roots.push_back(bisectionmethod(vprev,ve,h));
        }
        vprev=ve;
        errfprev=errfnow;
    }//this gives me all roots in roots ok its called allroots
    return roots;
}
//vectors. store things. how how how
void saveitall(std::vector<double> themroots, double stepsolution){//void bc its just to store does nothin
//i want a csv with all the roots to be able toplot solutions!! 
//how. solutions are all in rk whih is a struct
//store vector be careful, recall anything not in main this is general. in main i specify. BE CAREFUL!
std::vector<allsolns> savesolns;// EMPTY list to store the solutions, put in for let it update. 
//need all my solns
for (int k=0; k< themroots.size(); k++){
        //this just swipes all the roots we have and will save them. but i need the rk trajectory.
        savesolns.push_back(rungekutta(1,0,15,stepsolution,themroots[k]));//careful, h now is stepsoln, v0 is the roots so them roots. these parameters are a liiiitle confusing at this point.
    }//i pushed back now i need file similar to assignment 1, not too bad. but now i put outside loop bc i want to save 13 different roots. 
    //nevermind, i will make 13 csvs and put them in the same python graph. 
   
    for (int l=0; l< themroots.size(); l++){ //l goes over the 13 solns
        std:: string graphname = "solutionnumber_"+ std::to_string(l) +".csv";
        std::ofstream file(graphname); //saves in csv for the graph in python
         file<< "t,x,v\n";// this is to put it in ordee
         std::setprecision(10);
         for (int i=0; i<savesolns[l].tva.size(); i++) { // we create a loop that goes over the whole vectors
        file << savesolns[l].tva[i] << "," << savesolns[l].xva[i] << ","<< savesolns[l].vva[i]<< "\n"; } // everything with what corresponds.
        file.close(); //close. now what i goes over the stepsiz. :D
         }
}

//TRY TO DO CONVERGENCE ANALYSIS

int main(){
    double h=0.005;
    double root =bisectionmethod (1.0, 1.3, h);
    std::cout <<std::setprecision(8) <<"The slope x'(0) is "<< root<<"\n";
    std::cout<< "Can be verified because x(15) gives"<<"\n";

    double xfifteen= rungekutta (1.0, 0.0, 15.0,h,root).xva.back();//xva back fixes incompatibility bc it gives the 1 number we want not whole struct 3 things in my casee
    std::cout<<xfifteen<<"\n";
    
        // define roots vector again bc it was just in allroots thing
    std::vector<double> roots=allroots();
    std::cout <<"There are"<<roots.size()<<"roots"<<"\n";
    std::cout <<"Which are"<<"\n";
    for (int k=0; k< roots.size(); k++){
        std::cout<<"x'(0)"<<std::setprecision(10)<<roots[k]<< "\n";
    }
    saveitall(roots, h); //this simply generates them files for them roots from 0 to 12. its ok. 
    return 0;


}
