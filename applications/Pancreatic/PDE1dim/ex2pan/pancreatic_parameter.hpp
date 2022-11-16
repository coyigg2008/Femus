#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include "MultiLevelSolution.hpp"

using namespace femus;

// Table 2 Parameters descriptions and values
double D_D = 8.64e-7;
double D_C = 8.64e-7;
double D_T4 = 8.64e-7;
double D_T8 = 8.64e-7;
double D_Tr = 8.64e-7;
double D_Talpha = 1.29e-2;
double D_I2 = 9.95e-2;
double D_I12 = 6.04e-2;
double D_m1 = 0.13;
double D_m2 = 0.13;
double D_S = 9.26e-2;
double D_EC = 1.23e-4;
double D_ED = 1.23e-4;
double lambda_C = 0.46;
double lambda_1 = 3.12;
double lambda_T8C = 1.2e2;
double lambda_EC = 4.85e-8;
double lambda_ED = 5.23e-3;
double lambda_S = 1.9;
double lambda_D = 8.0e-7;
double lambda_T4 = 5;
double lambda_T4I2 = 0.25;
double lambda_T8 = 3.5;
double lambda_T8I2 = 0.25;
double lambda_Tr = 2.5;
double lambda_I12D = 5e-5;
double lambda_TalphaD = 1.44e-5;
double lambda_I2T4 = 3.0e-8;
double lambda_TalphaED = 0.123;
double lambda_m1EC = 2.22e-3;
double lambda_m2EC = 1.32;
double lambda_SED = 15.2;

//Table 3 Parameters descriptions and values
double d_T4 = 1.97e-1;
double d_T8 = 1.8e-1;
double d_Tr = 1.0e-1;
double d_C = 0.17;
double d_D = 0.1;
double d_Talpha = 35;
double d_I2 = 2.376;
double d_I12 = 1.38;
double d_EC = 21.8;
double d_ED = 21.8;
double d_m1 = 0.5;
double d_m2 = 0.5;
double d_S = 15;
double C_0 = 0.8;
double K_C = 0.4;
double K_T4 = 2e-3;
double K_T8 = 2e-3;
double K_Tr = 5e-4;
double T_0 = 1.2e-3;
double T_80 = 8e-4;
double T_hat4 = 3e-3;
double T_hat8 = 4e-3;
double K_I2 = 2.37e-11;
double K_I12 = 1.5e-10;
double K_m1 = 2.e-12;
double K_m2 = 1.2e-12;
double K_EC = 4.5e-10;
double K_D = 4e-6;
double K_S = 1e-10;
double K_Talpha = 1e-12;
double alpha = 1.0;
double A_S = 1.9e-9;



