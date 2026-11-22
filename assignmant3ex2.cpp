
#include <iostream> //general
#include <cstdlib>//for random numbers
#include <ctime> //for random numbers
#include <cmath> // math functions like exp, pow, etc
#include <random>// distributions etc
#include <algorithm> // for std::max, std::min, etc
#include <vector> // to use vectors, we need to include this library
#include <fstream> // for files
#include <iomanip> //for precision 

//laplace equation must not be defined because the SOR method contains it implicitly. the method is built around that equation

const int N=150; // we are defining dimensions, need to put const otherwise for the array phi it takes it as a variable and thats not possible
double h =1.0/N;
double phi [N+1][N+1]; //N+1 bc dividing in 100 sections gives 101 sections.
bool fixpoints [N+1][N+1]; // true is boundary or false is not. if false w appy the laplac


//separate and get the boundary points and the other points
//need a boundary values function but first we need to know which points are boundary points
//and which are not. 
//construct the grid of the image. first empty. then put the info inside

void grid(){

    for (int i=0; i<N+1;i++) { //phi sub i j. i is one dir j the other. go over all directions. two loops, one of i other of j
        for (int j=0;j<N+1; j++){
             phi [i][j]=0.0;
            fixpoints[i][j]= false;

        }
    }

}

void boundaries(){ //normal b
    //divide in i and j

    for (int i=0; i<N+1;i++){
        double x=i*h;
        phi[i][N] =x*x;    
        fixpoints[i][N]= true;
    }
    for (int j=0; j<N+1; j++){
        double y=h*j;
        phi[N][j] =y*y;    
        fixpoints[N][j]= true;
// this is for the fixed points, we delete the previous ones from pqrt 1 bc now we modify those
    }// neumann conditions are not fixed points, are a condition of the calculus. so the new condits. will go in the SOR phase
}

void extras(){ //box and line
    //box
    int boxx1=(int)(0.2 * N + 0.5);
    int boxx2= (int)(0.6 * N + 0.5);//0.5 to round up.
    int boxy1 =(int)(0.7 * N + 0.5);
    int boxy2=(int)(0.9 * N + 0.5);
// has to go over all these points
    for (int i = boxx1; i <= boxx2; i++) {
        for (int j = boxy1; j <= boxy2; j++) {
            phi[i][j] = 1.0;
            fixpoints[i][j] = true;
        }
    }

    //L Line horiz
    int linx1 = (int)(0.6 * N + 0.5);
    int linx2 = (int)(0.8 * N + 0.5);
    int height_linx  = (int)(0.5 * N + 0.5);
    //l line vertival
    int width_liny  = (int)(0.8 * N + 0.5);
    int liny1 = (int)(0.1 * N + 0.5);
    int liny2 = (int)(0.5 * N + 0.5);
    //go through all values
    for (int i = linx1; i <= linx2; i++) {
        phi[i][height_linx] = -1.0;
        fixpoints[i][height_linx] = true;

    }

    for (int j = liny1; j <= liny2; j++) {
        phi[width_liny][j] = -1.0;
        fixpoints[width_liny][j] = true;
    }


}
//now sor 
int sor(double w, double tolerance){ //now we start at 0 not 1
int sites=0; //n times we do loop
double maxdiffe =1.0; //itds ok as long as its bihher thean tolerancde
while (maxdiffe > tolerance) { // repeat until we get ehat we want
    maxdiffe=0.0; // will go updating
    for (int i = 0; i < N; i++) { //skip border points bc dont contribut
            for (int j = 0; j < N; j++) { //there is still a fixed border on the otherside so its N not N+1!!!!
                if (fixpoints [i][j]==true){ //want to compare not assifn
                    continue;
                } 
                //neuman
                double up;
                double right;
                double down;
                double left;

                right = phi[i+1][j];
                up = phi[i][j+1];

                if (i==0){
                    left= phi[1][j];}
                    else {
                        left= phi[i-1][j];
                    }
                if (j==0){
                    down=phi[i][1];

                }
                else{
                    down=phi[i][j-1];
                }


                double sorry = (1.0 - w) * phi[i][j]+ (w / 4.0) * (up+down+left+right);
                double change = fabs(sorry-phi [i][j]);
                if (change>maxdiffe){ // we wanr biggest change in chain, otherwise converfence wouldnt be correct
                    maxdiffe=change;
//sorry is sor but i had to change it to not have it twice 
                }
            phi[i][j] = sorry;
        } 
   
    }


    sites++; //update n of times we do the loop[]
}
return sites;
}
double derival() {
    int vali = (int)(0.3 * N + 0.5);  
    int valj = (int)(0.5 * N + 0.5);   

    double result = (phi[vali+1][valj] - phi[vali-1][valj]) / (2.0 * h); //divide by distn. we asked deriv we asjed slope so need to look at jjalf
    return result;
}

void graph(){
    std::ofstream file("phi2.csv");
    file << std::setprecision(10);   
    file << "x,y,phi\n";
    for (int i=0; i<N+1; i++){
        for (int j=0;j<N+1;j++){
            double x = i * h;
            double y = j * h;
            file << x << "," << y << "," << phi[i][j] << "\n";
        }
    }
    file.close();
}

double bestw(double tolerance){
    double wtrial = 1.0;
    double itemin = 1000000000.0;

    
    for (double w = 1.0; w < 1.99; w += 0.05){
        grid();
        boundaries();
        extras();
        int itenew = sor(w, tolerance);
        if (itenew < itemin){
            itemin = itenew;
            wtrial = w;
        }
    }

    
    double center = wtrial;
    for (double e = center - 0.05; e <= center + 0.05; e += 0.002){
        if (e < 1.0){
            continue;
        }
        if (e >= 1.999){
            break;
        }
        grid();
        boundaries();
        extras();
        int itenewsmol = sor(e, tolerance);
        if (itenewsmol < itemin){
            itemin = itenewsmol;
            wtrial = e;
        }
    }

    return wtrial;
}


int main(){
    double tol=0.00000000001;
    double yolo=bestw(tol);
    std::cout<<"The optimal w value is " <<yolo<<"\n";
    grid();
    boundaries();
    extras();
    int soror =sor(yolo,tol);
    graph();
std::cout << "Iterations = " << soror << "\n";
std::cout << "dphi/dx = " << derival() << "\n";
    
    


    return 0;

}