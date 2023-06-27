#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>
//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  V_skin - lx3 matrix of vertices of the display mesh
//Output:
//  N - the lxn sparse skinning matrix
bool IsInside(Eigen::Vector4d &phi,Eigen::Ref<const Eigen::MatrixXd> V, const Eigen::VectorXi& element,
              const Eigen::Vector3d& vertex)
{
    phi_linear_tetrahedron(phi,V,element,vertex);
    bool result=true;
    for(int i=0;i<4;++i)
    {
        if(phi[i]<0.0||phi[i]>1.0)
        {
            result=false;
            break;
        }
    }
    return result;
}

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {
    unsigned int n_tet=T.rows();
    unsigned int n_skin_verts=V_skin.rows();
    unsigned int n_tet_verts=V.rows();
    N.resize(n_skin_verts,n_tet_verts);
    N.setZero();
    std::vector<Eigen::Triplet<double>> TripletList;
    for(int j=0;j<n_skin_verts;++j)
    {

        Eigen::Vector3d vertex=V_skin.row(j);
        for(int i=0;i<n_tet;++i)
        {
            Eigen::VectorXi element=T.row(i);
            Eigen::Vector4d phi;

            if(IsInside(phi,V,element,vertex))
            {
                for(int k=0;k<4;++k)
                {
                    TripletList.emplace_back(j,element[k],phi[k]);
                }
                break;
            }
        }
    }

    N.setFromTriplets(TripletList.begin(),TripletList.end());
    N.makeCompressed();

    
}