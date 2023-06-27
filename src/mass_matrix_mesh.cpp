#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>
//Input:
//  qdot - generalized velocity for the FEM system
//  T - the mx4 vertex indices for tet mesh
//  density - density of material
//  v0 - the undeformed tetrahedra volumes
//Output:
//  M - Sparse mass matrix for the whole mesh.
void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {
    unsigned int n_tet=T.rows();
    unsigned int q_size=qdot.size();
    M.resize(q_size,q_size);
    M.setZero();
    for(int n=0;n<n_tet;++n)
    {
        unsigned int ids[4];
        Eigen::RowVectorXi element(4);
        for(int i=0;i<4;++i)
        {
            ids[i]=T(n,i);
        }
        element=T.row(n);
        Eigen::SparseMatrixd P(12,q_size);
        P.setZero();
        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.clear();
        for(int i=0;i<4;++i)
        {
            for(int j=0;j<3;++j)
            {
                tripletList.emplace_back(3*i+j,3*element(i)+j,1.0);
            }
        }
        P.setFromTriplets(tripletList.begin(),tripletList.end());
        P.makeCompressed();
        Eigen::Matrix1212d Me;
        mass_matrix_linear_tetrahedron(Me,qdot,element,density,v0[n]);
        M+=(P.transpose()*Me*P).sparseView();

    }
    

}