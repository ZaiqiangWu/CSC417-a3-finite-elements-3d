#include <assemble_stiffness.h>
#include <iostream>
#include "timer.h"
//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  K - the sparse, global stiffness matrix
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D) {
    //std::cout<<"assemble_stiffness start"<<std::endl;
    unsigned int n_tet=T.rows();
    unsigned int q_size=qdot.size();
    K.resize(q_size,q_size);
    K.setZero();
    std::vector<Eigen::Triplet<double>> TripletList;
    Timer timer=Timer();
    timer.begin();

    for(int n=0;n<n_tet;++n)
    {

        Eigen::RowVectorXi element(4);
        element=T.row(n);
        Eigen::Matrix1212d Ke;

        d2V_linear_tetrahedron_dq2(Ke,q,V,element,v0[n],C,D);
        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.clear();
        for(int i=0;i<4;++i)
        {
            for(int j=0;j<4;++j)
            {
                for(int a=0;a<3;++a)
                {
                    for(int b=0;b<3;++b)
                    {
                        TripletList.emplace_back(3*element(i)+a,3*element(j)+b,-Ke(3*i+a,3*j+b));
                    }

                }
                //tripletList.emplace_back(3*i+j,3*element(i)+j,1.0);
            }
        }

    }

    std::cout<<"H for single item:"<<timer.get_elapsed_time()<<std::endl;

    K.setFromTriplets(TripletList.begin(),TripletList.end());
    K.makeCompressed();
    //std::cout<<"H for single item:"<<timer.get_elapsed_time()<<std::endl;

    //std::cout<<"assemble_stiffness end"<<std::endl;
    };
