#include <fstream>
#include <iostream>


template <typename Vect, typename ParMat, typename ParVect>
class melt_sphere {
    private:
        int M; // number of space steps (-)
        int N; // number of time steps (-)
        double dx; // non-dimensional space step (-)
        double dt; // non-dimensional time step (-)
        double dtheta; // non-dimensional semi-interval for mushy ice (-)
        double end_time; // end of simulation 
        double Te; // non-dimensional external temperature (-) 
        double stefan_number; // Stefan number (-)
        ParVect param_beta; // array of normalized parameter beta (-)
        Vect a,b,c,f; // three-diagonal and right-end vectors
        // private methods
        void generate_sparse_matrix(int ts);
        void generate_right_hand(int ts);
        void generate_beta(int ts);
        double beta(double theta); // help function for generate_beta
        void lu_factorization(); 
        void backward();
        void forward();
    
    public:
        melt_sphere(int m, int n, 
                         double et, double deltaT, double stef);
        ParMat sol; // matrix that store the temperature solution
        // note that this is for illustration purpose only because
        // it necessitates a lot of memory! DO NOT do such things in
        // production 
        void save_in_file(std::string file_name);
        void solve();
};

template <typename Vect, typename ParMat, typename ParVect>
melt_sphere<Vect,ParMat,ParVect>::melt_sphere(int m, int n, 
                         double et, double deltaT, double stef){
    M = m;
    N = n;
    dx = 1./m;
    end_time = et;
    dt = end_time/n;
    dtheta = deltaT;
    Te = 0.;
    ParMat solTemp(n+1, Vect(m+1));
    sol = solTemp;
    // HACK: here we add a 1.01 to dtheta in order to make sure the 
    // initial state is solid and not mushy
    std::fill(sol[0].begin(), sol[0].end(), 1.+dtheta*1.01);
    Vect aTemp(M,1.);
    a=aTemp;
    Vect bTemp(M+1,1.);
    b=bTemp;
    Vect cTemp(M,1.);
    c=cTemp;
    Vect fTemp(M+1,0.);
    f=fTemp;
    param_beta.resize(M+1);
    stefan_number = stef;
}

template <typename Vect, typename ParMat, typename ParVect>
void melt_sphere<Vect,ParMat,ParVect>::save_in_file(std::string file_name){
      std::ofstream temp_file(file_name + "_temp.txt", std::ios::out | std::ios::trunc);  
      if(temp_file)  
      { 
          for(int i=0; i<M+1; i++){
              for(int j=0; j<N+1; j++){
                  temp_file << " " << sol[j][i];
              }
              temp_file << std::endl;
          }
        temp_file.close();
        }
      else 
      {
        std::cerr << "Error while opening the text file." << std::endl;
      }
}

template <typename Vect, typename ParMat, typename ParVect>
void melt_sphere<Vect,ParMat,ParVect>::generate_sparse_matrix(int ts){
    ParVect r = param_beta*dt/(2.*dx*dx);
    b[0] = 1.0;
    c[0] = -1.0;
    for(int i=1; i<M; i++){
        a[i-1] = -r[i]*(1.-1./i);
        b[i] = 1.+2.*r[i];
        c[i] = -r[i]*(1.+1./i);
    }
    a[M-1]=0.;
    b[M] =1.;
    return;
}

template <typename Vect, typename ParMat, typename ParVect>
void melt_sphere<Vect,ParMat,ParVect>::generate_right_hand(int ts){
    ParVect r = param_beta*dt/(2.*dx*dx);
    f[0]=0.0;
    f[1]=(1.-2.*r[1])*sol[ts][1]+r[1]*(1.+1./1.)*sol[ts][2];
    for(int i=2; i<M; i++){
        f[i]=r[i]*(1.-1./i)*sol[ts][i-1]+(1.-2.*r[i])*sol[ts][i]+r[i]*(1.+1./i)*sol[ts][i+1];
    }
    f[M] = Te;
    return;
}

template <typename Vect, typename ParMat, typename ParVect>
double melt_sphere<Vect,ParMat,ParVect>::beta(double theta){
    double res;
    if(theta > 1.+dtheta){
        res = 1.;
    }
    else if(theta < 1.-dtheta){
        res = 1.;
    }
    else{
        res = 1./(1.+stefan_number/(2.*dtheta));
    }
    return res;
}

template <typename Vect, typename ParMat, typename ParVect>
void melt_sphere<Vect,ParMat,ParVect>::generate_beta(int ts){
    for(int i=0; i<M; i++){
        param_beta[i] = melt_sphere<Vect,ParMat,ParVect>::beta(sol[ts][i]);
    }
    return;
}

template <typename Vect, typename ParMat, typename ParVect>
void melt_sphere<Vect,ParMat,ParVect>::lu_factorization(){
    c[0]=c[0]/b[0];
    b[0]=b[0];
    for(int i=1; i<M; i++){
        b[i]=b[i]-a[i-1]*c[i-1];
        c[i]=c[i]/b[i];
    }
    b[M]=b[M]-a[M-1]*c[M-1];
}

template <typename Vect, typename ParMat, typename ParVect>
void melt_sphere<Vect,ParMat,ParVect>::forward(){
    f[0]=f[0]/b[0];
    for(int i=1; i<M+1; i++){
        f[i]=(f[i]-a[i-1]*f[i-1])/b[i];
    }
}

template <typename Vect, typename ParMat, typename ParVect>
void melt_sphere<Vect,ParMat,ParVect>::backward(){
    for(int i=M-1; i>=0; i--){
        f[i]=f[i]-c[i]*f[i+1];
    }
}

template <typename Vect, typename ParMat, typename ParVect>
void melt_sphere<Vect,ParMat,ParVect>::solve(){
    for(int i=1; i<=N; i++){
        melt_sphere<Vect,ParMat,ParVect>::generate_beta(i-1);
        melt_sphere<Vect,ParMat,ParVect>::generate_right_hand(i-1);
        melt_sphere<Vect,ParMat,ParVect>::generate_sparse_matrix(i-1);
        melt_sphere<Vect,ParMat,ParVect>::lu_factorization();
        melt_sphere<Vect,ParMat,ParVect>::forward();
        melt_sphere<Vect,ParMat,ParVect>::backward();
        sol[i]=f;
    }
}
