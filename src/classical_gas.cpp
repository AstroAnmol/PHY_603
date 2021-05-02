#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <random>
#include "classical_gas.h"
//double random number generator
double dRand(double dMin, double dMax){
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_real_distribution<double>    distr(dMin, dMax);
    return distr(generator);
}
// int random number generator
int iRand(int range_from, int range_to) {
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>    distr(range_from, range_to);
    return distr(generator);
}
// Intiator
classical_gas::classical_gas(int number, int iterations, double dens){
    n=number;
    iter=iterations;
    density=dens;
    V=n/density;
    L=std::pow(V,1.0/3.0);
    // std::cout<<L<<std::endl;
    delta=L/20.0;
    T=beta/kB;
    Virial_all=Eigen::ArrayXd::Zero(iter+1);
    Acceptance=Eigen::ArrayXd::Zero(iter);
    initiate_particles();
    Iter_Conf=Initial_Conf;
    Virial_all(0)=Virial_calc(Initial_Conf); //initial virial
    for (int i = 0; i < iter; i++){
        // select one particle randomly
        int sel=iRand(0, n-1);
        // calculate the configurtion energy
        double Iter_energy;
        Iter_energy =part_energy(Iter_Conf, sel);
        // Calculate proposal new position
        Eigen::Vector3d new_pos;
        new_pos=proposal_particle_disp(Iter_Conf.col(sel));
        // std::cout<<Iter_Conf.col(sel).transpose()<<'|'<<new_pos.transpose()<<std::endl;
        // Generate new configuration
        Eigen::ArrayXXd Propose_Conf;
        Propose_Conf=Iter_Conf;
        Propose_Conf.col(sel)=new_pos;
        double Propose_energy;
        Propose_energy=part_energy(Propose_Conf, sel);
        //std::cout<<Propose_energy<<' ';
        //std::cout<<Iter_energy<<' ';
        // Check acceptance
        if (std::exp(-beta*(Propose_energy-Iter_energy))>dRand(0,1)){
            Iter_Conf=Propose_Conf;
            Acceptance(i)=1;
        }
        //std::cout<<Acceptance(i)<<std::endl;
        Virial_all(i+1)=Virial_calc(Iter_Conf);
    }
    Equilibirum_Conf=Iter_Conf;
    std::cout<<"Acceptance Rate: "<<Acceptance.sum()/iter<<std::endl;
    Jackknife();
}

//relative distance between two particles
Eigen::Vector3d classical_gas::rel_dist(Eigen::Vector3d Ra, Eigen::Vector3d Rb){
    Eigen::Vector3d towardsa;
    for(int i=0; i<3; i++){
        double a=Ra(i), b=Rb(i);
        if(abs(a-b)<L/2.0){
            towardsa(i)=a-b;}
        else if(b>=a){
            towardsa(i)=a-b+L;}
        else{
            towardsa(i)=a-b-L;}
    }
    return towardsa;
}
// force on a due to b
Eigen::Vector3d classical_gas::force_gaussian(Eigen::Vector3d Ra, Eigen::Vector3d Rb){ 
    Eigen::Vector3d del_r = rel_dist(Ra, Rb);
    double rabs2 = del_r.dot(del_r);
    return del_r*(2.0*exp(-rabs2/R2))/R2;
}
// Potential between a and b
double classical_gas::potential_gaussian(Eigen::Vector3d Ra, Eigen::Vector3d Rb){ 
    Eigen::Vector3d del_r = rel_dist(Ra, Rb);
    double rabs2 = del_r.dot(del_r);
    return exp(-rabs2/R2);
}
// initiate particles randomly
void classical_gas::initiate_particles(){
    Eigen::ArrayXXd Init_Conf(3,n);
    for (int i = 0; i < n; i++){ // loop over n particles
        Eigen::Vector3d part_pos;
        for (int a=0; a<3; a++){ // loop over three position components
            part_pos(a)=dRand(0,L);
        }
        Init_Conf.col(i)=part_pos;
    }
    Initial_Conf=Init_Conf;
}
// calculate conf energy
double classical_gas::conf_energy(Eigen::ArrayXXd Conf){
    double energy=0;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i<j){
                energy=energy+potential_gaussian(Conf.col(i),Conf.col(j));
            }
        }
    }
    return energy;
}
// calculate energy of kth particle in a conf
double classical_gas::part_energy(Eigen::ArrayXXd Conf, int k){
    double energy=0;
    for (int i = 0; i < n; i++){
        if (i!=k){
            energy=energy+potential_gaussian(Conf.col(i),Conf.col(k));
        }
    }
    return energy;
}

// particle displacement proposal
Eigen::Vector3d classical_gas::proposal_particle_disp(Eigen::Vector3d org_pos){
    Eigen::Vector3d new_pos;
    for (int i = 0; i < 3; i++){
        new_pos(i)=delta*(dRand(0,1)-0.5);
    }
    return new_pos;
} 
// Calculate Virial
double classical_gas::Virial_calc(Eigen::ArrayXXd Conf){
    double vir=0;
    Eigen::Vector3d force;
    Eigen::Vector3d rel_distance;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i!=j){
                force=force_gaussian(Conf.col(i),Conf.col(j));
                rel_distance=rel_dist(Conf.col(i),Conf.col(j));
                vir=vir + force.dot(rel_distance);
            }
        }
    }
    vir=vir/2.0;
    return vir;
}

// Jackknife Analysis
void classical_gas::Jackknife(){
    //compute average exluding ith measurement
    Eigen::ArrayXd Avg_exd_i(iter+1);
    for (int i = 0; i < iter+1; i++){
        double avg;
        avg=Virial_all.sum()-Virial_all(i);
        avg=avg/(iter);
        Avg_exd_i(i)=avg;
    }
    // compute best estimate of virial
    Virial_avg=Avg_exd_i.sum()/(iter+1);
    // compute error bar
    double a;
    for (int i = 0; i < iter+1; i++){
        a = a + std::pow((Avg_exd_i(i) - Virial_avg),2);
    }
    error_bar=std::sqrt( iter*a/(iter+1) );
}   

//Get functions
double classical_gas::get_avg_Virial(){
    return Virial_avg;
}
double classical_gas::get_error_bar(){
    return error_bar;
}

void classical_gas::get_all_Virial(std::string name){
    Eigen::IOFormat csv(Eigen::FullPrecision, 0, ", ", "\n", "", "", "", "");
    std::ofstream theFile;
    std::ostringstream oss;
    oss <<name<<"_Virial_"<<n<<'_'<<density<<'_'<<beta<<'_'<<iter;
    std::string file_name;
    file_name=oss.str();
    theFile.open(file_name);
    theFile << "Virials" << std::endl;
    theFile << Virial_all.format(csv)<< std::endl;
    theFile.close();
}