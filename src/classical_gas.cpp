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
classical_gas::classical_gas(int number, int iterations){
    n=number;
    iter=iterations;
    V=n/density;
    L=std::pow(V,1.0/3.0);
    delta=L/20.0;
    T=beta/kB;
    initiate_particles();
    Iter_Conf=Initial_Conf;
    for (int i = 0; i < iter; i++){
        // select one particle randomly
        int sel=iRand(0, n-1);
        // calculate the configurtion energy
        double Iter_energy;
        Iter_energy =conf_energy(Iter_Conf);
        // Calculate proposal new position
        Eigen::Vector3d new_pos;
        new_pos=proposal_particle_disp(Iter_Conf.col(sel));
        // Generate new configuration
        Eigen::ArrayXXd Propose_Conf;
        Propose_Conf=Iter_Conf;
        Propose_Conf.col(sel)=new_pos;
        double Propose_energy;
        Propose_energy=conf_energy(Propose_Conf);
        // Check acceptance
        if (std::exp(-beta*(Propose_energy+Iter_energy)) >dRand(0,1)){
            Iter_Conf=Propose_Conf;
        }
    }
    Equilibirum_Conf=Iter_Conf;
    Virial_calc();
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
    double energy;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i!=j){
                energy=energy+potential_gaussian(Conf.col(i),Conf.col(j));
            }
        }
    }
    return energy/2.0;
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
void classical_gas::Virial_calc(){
    double vir;
    Eigen::Vector3d force;
    Eigen::Vector3d rel_distance;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i!=j){
                force=force_gaussian(Equilibirum_Conf.col(i),Equilibirum_Conf.col(j));
                rel_distance=rel_dist(Equilibirum_Conf.col(i),Equilibirum_Conf.col(j));
                vir=vir + force.dot(rel_distance);
            }
        }
    }
    vir=vir/2.0;
    Virial=vir;
}

//Get functions
double classical_gas::get_Virial(){
    return Virial;
}