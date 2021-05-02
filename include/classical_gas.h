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
        double density; // n/V
        double T; //Temperature
        double kB=1.380649E-23; // Boltzmann Constant
        double beta=1; // Boltzmann factor
        double delta; //delta
        Eigen::ArrayXd Virial_all; // Array of Virial for each configuration
        double Virial_avg; //average virial
        double error_bar; //error bar
        Eigen::ArrayXd Acceptance;
        //functions
        Eigen::Vector3d rel_dist(Eigen::Vector3d Ra, Eigen::Vector3d Rb);
        Eigen::Vector3d force_gaussian(Eigen::Vector3d Ra, Eigen::Vector3d Rb);
        double potential_gaussian(Eigen::Vector3d Ra, Eigen::Vector3d Rb); 
        void initiate_particles();
        double conf_energy(Eigen::ArrayXXd Conf);
        double part_energy(Eigen::ArrayXXd Conf, int k);
        Eigen::Vector3d proposal_particle_disp(Eigen::Vector3d org_pos);
        double Virial_calc(Eigen::ArrayXXd Conf);
        void Jackknife();
    public:
        classical_gas(int number, int iterations, double dens);
        Eigen::ArrayXXd Initial_Conf;
        Eigen::ArrayXXd Iter_Conf;
        Eigen::ArrayXXd Equilibirum_Conf;
        double get_avg_Virial();
        double get_error_bar();
        void get_all_Virial(std::string name); //makes a csv file for virials
};
#endif