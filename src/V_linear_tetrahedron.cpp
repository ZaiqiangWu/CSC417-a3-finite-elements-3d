#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X){
        Eigen::Matrix3d F,Dm,Cm;
        for(int i=0;i<3;++i)
        {
            Dm.col(i)=V.row(element(i+1))-V.row(element(0));
            Cm.col(i)=q.block(3*element(i+1),0,3,1)-q.block(3*element(0),0,3,1);
        }
        F=Cm*Dm.inverse();

        psi_neo_hookean(e,F,C,D);

       
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);  
    
}