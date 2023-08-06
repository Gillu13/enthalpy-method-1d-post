#include <vector>
#include <valarray>
#include "sphere_phase_change.cpp"

typedef std::vector<double> Vect;
typedef std::vector<std::vector<double>> ParMat;
typedef std::valarray<double> ParVect;

int main(int argc, char* argv[]){
int M, N;
double stef, k, et;

if(argc==6){
    M = std::stoi(argv[1]); // number of length steps (-)
    N = std::stoi(argv[2]); // number of time steps (-)
    // environmental conditions
    stef = std::stod(argv[3]);
    k = std::stod(argv[4]); // dtheta (-)
    et=std::stod(argv[5]); // end time (-)
}
else{
    M=400; // number of space steps (-)
    // environmental conditions
    stef = 1.;
    k = 0.1; // dtheta (-)
    // time information
    et=0.3; // end time (-)
    N=400; // number of time steps (-)
}

double end_time=et*stef; // end time (-)
double dTemp=k; // non-dimensional dtheta [-]

std::string file_name="solidification";


melt_sphere<Vect,ParMat,ParVect> prob(M, N, end_time, dTemp, stef);
prob.solve();
prob.save_in_file(file_name);

return 0;
}
