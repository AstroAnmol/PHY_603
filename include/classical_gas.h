#ifndef CLASSICAL_H
#define CLASSICAL_H
#include <eigen-3.3.7/Eigen/Dense>

class classical_gas{
    private:
        //constants
        int n; // no of particles
        int iter; // no of iterations
        double L; // Length of the box (m)
        double V; // Volume of the box (m^3)  
        double R2=1.0; // constant for the gaussian potential
        double density=0.16; // n/V
        double T; //Temperature
        double kB=1.380649E-23; // Boltzmann Constant
        double beta=1; // Boltzmann factor
        double delta; //delta
        double Virial; //virial
        //functions
        Eigen::Vector3d rel_dist(Eigen::Vector3d Ra, Eigen::Vector3d Rb);
        Eigen::Vector3d force_gaussian(Eigen::Vector3d Ra, Eigen::Vector3d Rb);
        double potential_gaussian(Eigen::Vector3d Ra, Eigen::Vector3d Rb); 
        void initiate_particles();
        double conf_energy(Eigen::ArrayXXd Conf);
        Eigen::Vector3d proposal_particle_disp(Eigen::Vector3d org_pos);
        void Virial_calc();
    public:
        classical_gas(int number, int iterations);
        Eigen::ArrayXXd Initial_Conf;
        Eigen::ArrayXXd Iter_Conf;
        Eigen::ArrayXXd Equilibirum_Conf; 
        double get_Virial();
};
#endif