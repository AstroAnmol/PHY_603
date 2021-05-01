#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <fstream>
#include <eigen-3.3.7/Eigen/Dense>
#include <cmath>
#include <string>
#include <random>
#include "classical_gas.h"

int main(){
    int n=40;
    int iterations=10000;
    classical_gas gas(n, iterations);

    double Virial;
    Virial = gas.get_Virial();
    std::cout<<Virial/n<<std::endl; 
}
