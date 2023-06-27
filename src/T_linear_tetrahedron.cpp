#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
    Eigen::Matrix43d dphi;
    Eigen::Matrix1212d I=Eigen::Matrix1212d::Identity();
    Eigen::Matrix1212d Me;
    mass_matrix_linear_tetrahedron(Me,qdot,element,density,volume);
    T=0.5*qdot.transpose()*Me*qdot;
}