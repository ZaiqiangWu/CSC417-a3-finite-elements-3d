#include <assemble_forces.h>
#include <iostream>
//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system
void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) {

    unsigned int n_tet=T.rows();
    unsigned int q_size=q.size();
    f.resize(q_size);
    f.setZero();

    for(int n=0;n<n_tet;++n)
    {

        Eigen::RowVectorXi element(4);

        element=T.row(n);
        Eigen::Vector12d fe;
        //std::cout<<q<<std::endl;
        dV_linear_tetrahedron_dq(fe,q,V,element,v0[n],C,D);
        fe=-fe;
        //std::cout<<fe.hasNaN()<<std::endl;

        for(int i=0;i<4;++i)
        {
            for(int j=0;j<3;++j)
            {
                f(3*element(i)+j)+=fe(3*i+j);
            }
        }
    }

    };