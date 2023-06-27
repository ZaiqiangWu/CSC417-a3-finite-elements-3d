#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void compute_dFdq(Eigen::MatrixXd &dFdq,Eigen::Ref<const Eigen::Matrix3d> inv_Dm)
{
    Eigen::Matrix43d assemble=Eigen::Matrix43d::Ones();
    assemble*=-1.0;
    assemble.block(1,0,3,3)=Eigen::Matrix3d::Identity();
    Eigen::Matrix43d C=assemble*inv_Dm;
    dFdq.resize(9,12);
    for(int i=0;i<3;++i)
    {
        for(int j=0;j<3;++j)
        {
            for(int m=0;m<3;++m)
            {
                for(int n=0;n<4;++n)
                {
                    if(i!=m)
                    dFdq(3*j+i,3*n+m)=0.0;
                    else
                        dFdq(3*j+i,3*n+m)=C(n,j);
                }

            }
        }
    }
}

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
       Eigen::Matrix3d F,Dm,Cm;
       for(int i=0;i<3;++i)
       {
           Dm.col(i)=V.row(element(i+1))-V.row(element(0));
           Cm.col(i)=q.block(3*element(i+1),0,3,1)-q.block(3*element(0),0,3,1);
       }
       Eigen::Matrix3d inv_Dm=Dm.inverse();
       F=Cm*inv_Dm;
       //lot of zero in q??!!
       //std::cout<<q<<std::endl;
       //std::cout<<"a"<<std::endl;

       Eigen::Vector9d dpsidF;
       dpsi_neo_hookean_dF(dpsidF,F,C,D);

       Eigen::MatrixXd dFdq;
       compute_dFdq(dFdq,inv_Dm);
       dV=dFdq.transpose()*dpsidF;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);


}