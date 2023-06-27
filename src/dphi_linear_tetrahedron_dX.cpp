#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>
using namespace Eigen;
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    Vector3d x[4];
    for(int i=0;i<4;++i)
        x[i]=V.row(element(i));
    for(int i=0;i<4;++i)
    {
        Vector3d N=(x[(i+1)%4]-x[(i+2)%4]).cross(x[(i+1)%4]-x[(i+3)%4]).normalized();
        Vector3d a=N/((x[i+1]-x[i]).dot(N));
        //phi[i]=(x-X[i]).dot(N)/((X[i+1]-X[i]).dot(N));
        dphi.row(i)=a;
    }
}