#include <dpsi_neo_hookean_dF.h>
#include <iostream>
//Input:
//  F - the dense 3x3 deformation gradient
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  psi - the 9x1 gradient of the potential energy wrt to the deformation gradient
void dpsi_neo_hookean_dF(Eigen::Vector9d &dw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
    double f00=F(0,0);
    double f01=F(0,1);
    double f02=F(0,2);
    double f10=F(1,0);
    double f11=F(1,1);
    double f12=F(1,2);
    double f20=F(2,0);
    double f21=F(2,1);
    double f22=F(2,2);
    //Colum major
    dw(0)= C*(2*f00/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 2.0/3.0) + (-2.0/3.0*f11*f22 + (2.0/3.0)*f12*f21)*(pow(f00, 2) + pow(f01, 2) + pow(f02, 2) + pow(f10, 2) + pow(f11, 2) + pow(f12, 2) + pow(f20, 2) + pow(f21, 2) + pow(f22, 2))/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 5.0/3.0)) + D*(2*f11*f22 - 2*f12*f21)*(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20 - 1) ;
    dw(1)= C*(2*f10/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 2.0/3.0) + ((2.0/3.0)*f01*f22 - 2.0/3.0*f02*f21)*(pow(f00, 2) + pow(f01, 2) + pow(f02, 2) + pow(f10, 2) + pow(f11, 2) + pow(f12, 2) + pow(f20, 2) + pow(f21, 2) + pow(f22, 2))/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 5.0/3.0)) + D*(-2*f01*f22 + 2*f02*f21)*(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20 - 1) ;
    dw(2)= C*(2*f20/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 2.0/3.0) + (-2.0/3.0*f01*f12 + (2.0/3.0)*f02*f11)*(pow(f00, 2) + pow(f01, 2) + pow(f02, 2) + pow(f10, 2) + pow(f11, 2) + pow(f12, 2) + pow(f20, 2) + pow(f21, 2) + pow(f22, 2))/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 5.0/3.0)) + D*(2*f01*f12 - 2*f02*f11)*(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20 - 1) ;
    dw(3)= C*(2*f01/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 2.0/3.0) + ((2.0/3.0)*f10*f22 - 2.0/3.0*f12*f20)*(pow(f00, 2) + pow(f01, 2) + pow(f02, 2) + pow(f10, 2) + pow(f11, 2) + pow(f12, 2) + pow(f20, 2) + pow(f21, 2) + pow(f22, 2))/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 5.0/3.0)) + D*(-2*f10*f22 + 2*f12*f20)*(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20 - 1) ;
    dw(4)= C*(2*f11/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 2.0/3.0) + (-2.0/3.0*f00*f22 + (2.0/3.0)*f02*f20)*(pow(f00, 2) + pow(f01, 2) + pow(f02, 2) + pow(f10, 2) + pow(f11, 2) + pow(f12, 2) + pow(f20, 2) + pow(f21, 2) + pow(f22, 2))/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 5.0/3.0)) + D*(2*f00*f22 - 2*f02*f20)*(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20 - 1) ;
    dw(5)= C*(2*f21/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 2.0/3.0) + ((2.0/3.0)*f00*f12 - 2.0/3.0*f02*f10)*(pow(f00, 2) + pow(f01, 2) + pow(f02, 2) + pow(f10, 2) + pow(f11, 2) + pow(f12, 2) + pow(f20, 2) + pow(f21, 2) + pow(f22, 2))/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 5.0/3.0)) + D*(-2*f00*f12 + 2*f02*f10)*(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20 - 1) ;
    dw(6)= C*(2*f02/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 2.0/3.0) + (-2.0/3.0*f10*f21 + (2.0/3.0)*f11*f20)*(pow(f00, 2) + pow(f01, 2) + pow(f02, 2) + pow(f10, 2) + pow(f11, 2) + pow(f12, 2) + pow(f20, 2) + pow(f21, 2) + pow(f22, 2))/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 5.0/3.0)) + D*(2*f10*f21 - 2*f11*f20)*(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20 - 1) ;
    dw(7)= C*(2*f12/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 2.0/3.0) + ((2.0/3.0)*f00*f21 - 2.0/3.0*f01*f20)*(pow(f00, 2) + pow(f01, 2) + pow(f02, 2) + pow(f10, 2) + pow(f11, 2) + pow(f12, 2) + pow(f20, 2) + pow(f21, 2) + pow(f22, 2))/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 5.0/3.0)) + D*(-2*f00*f21 + 2*f01*f20)*(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20 - 1) ;
    dw(8)= C*(2*f22/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 2.0/3.0) + (-2.0/3.0*f00*f11 + (2.0/3.0)*f01*f10)*(pow(f00, 2) + pow(f01, 2) + pow(f02, 2) + pow(f10, 2) + pow(f11, 2) + pow(f12, 2) + pow(f20, 2) + pow(f21, 2) + pow(f22, 2))/pow(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20, 5.0/3.0)) + D*(2*f00*f11 - 2*f01*f10)*(f00*f11*f22 - f00*f12*f21 - f01*f10*f22 + f01*f12*f20 + f02*f10*f21 - f02*f11*f20 - 1) ;

    //std::cout<<"F:"<<F.hasNaN()<<std::endl;
    //std::cout<<F.maxCoeff()<<std::endl;





}