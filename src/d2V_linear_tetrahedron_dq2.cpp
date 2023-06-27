#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>
#include "dV_linear_tetrahedron_dq.h"
#include "timer.h"

void d2V_linear_tetrahedron_dq2(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Matrix1212d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        
      //Code to compute non-integrated hessian matrix goes here
       Eigen::Matrix3d F,Dm,Cm;
       for(int i=0;i<3;++i)
       {
           Dm.col(i)=V.row(element(i+1))-V.row(element(0));
           Cm.col(i)=q.block(3*element(i+1),0,3,1)-q.block(3*element(0),0,3,1);
       }
       Eigen::Matrix3d inv_Dm=Dm.inverse();
       F=Cm*inv_Dm;
       Eigen::Matrix99d d2psidF2;
       d2psi_neo_hookean_dF2(d2psidF2,F,C,D);
       Eigen::MatrixXd dFdq;
       compute_dFdq(dFdq,inv_Dm);
       Eigen::Vector12d dFdq_vec;
       for(int r=0;r<4;++r)
       {
           for(int c=0;c<3;++c)
           {
               dFdq_vec[4*r+c]=dFdq(r,c);
           }
       }
       dV=dFdq.transpose()*d2psidF2*dFdq;
     
    };

    //integrate the non-integrated hessian across the tetrahedral element
    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);  
    

    //DO NOT REMOVE THIS CODE This code ensures that the hessian matrix is symmetric postive definite by projecting all
    //negative eigenvalues to small, postive values.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 12; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();

}
