#include <phi_linear_tetrahedron.h>
using namespace Eigen;
void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {
    Vector3d X[4];
    for(int i=0;i<4;++i)
        X[i]=V.row(element(i));
    for(int i=0;i<4;++i)
    {
        Vector3d N=(X[(i+1)%4]-X[(i+2)%4]).cross(X[(i+1)%4]-X[(i+3)%4]).normalized();
        phi[i]=(x-X[(i+1)%4]).dot(N)/((X[i]-X[(i+1)%4]).dot(N));
    }
}