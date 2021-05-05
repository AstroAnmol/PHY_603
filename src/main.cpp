#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <fstream>
#include <eigen-3.3.7/Eigen/Dense>
#include <cmath>
#include <string>
#include <random>
#include "classical_gas.h"


// //relative distance between two particles
// Eigen::Vector3d rel_dist(Eigen::Vector3d Ra, Eigen::Vector3d Rb){
//     Eigen::Vector3d towardsa;
//     for(int i=0; i<3; i++){
//         double a=Ra(i), b=Rb(i);
//         if(abs(a-b)<10/2.0){
//             towardsa(i)=a-b;}
//         else if(b>=a){
//             towardsa(i)=a-b+10;}
//         else{
//             towardsa(i)=a-b-10;}
//     }
//     return towardsa;
// }
// // force on a due to b
// Eigen::Vector3d force_gaussian(Eigen::Vector3d Ra, Eigen::Vector3d Rb){ 
//     Eigen::Vector3d del_r = rel_dist(Ra, Rb);
//     double rabs2 = del_r.dot(del_r);
//     return del_r*(2.0*exp(-rabs2/1))/1;
// }
// // Potential between a and b
// double potential_gaussian(Eigen::Vector3d Ra, Eigen::Vector3d Rb){ 
//     Eigen::Vector3d del_r = rel_dist(Ra, Rb);
//     double rabs2 = del_r.dot(del_r);
//     return exp(-rabs2/1);
// }

int main(){   
    // Eigen::Vector3d Ra, Rb;
    // Ra<<1,2,3;
    // Rb<<4,5,6;
    // std::cout<<potential_gaussian(Ra,Rb)<<std::endl;
    // std::cout<<potential_gaussian(Rb,Ra)<<std::endl;
    // std::cout<<force_gaussian(Ra,Rb)<<std::endl;
    // std::cout<<force_gaussian(Rb,Ra)<<std::endl;
    Eigen::ArrayXd n(17);
    n<<40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200;
    int iterations=10000;
    //Eigen::ArrayXd density(9);
    //density<<0.001,0.005,0.01,0.05,0.2,0.5,1,10,100;
    double density=0.1;
    //for (int j = 0; j < 9; j++){
        std::cout<<"Density: "<<density<<std::endl;
        std::cout<<"Beta: 1"<<std::endl;
        std::cout<<"Iterations: 10000"<<std::endl;
        for (int i = 0; i < 17; i++){
            std::cout<<"---------------------------------"<<std::endl;
            std::cout<<"Number of Particles: "<< n(i)<<std::endl;
            std::cout<<"---------------------------------"<<std::endl;
            classical_gas gas(n(i), iterations,density);
            double Virial, error;
            Virial = gas.get_avg_Virial();
            error = gas.get_error_bar();
            std::cout<<"Virial per particle: "<<Virial/n(i)<<std::endl;
            std::cout<<"Error bar: "<<error/n(i)<<std::endl<<std::endl;
            gas.get_all_Virial("Trial1");
        }
    //}
}