 
 #include <mass_matrix_linear_tetrahedron.h>

 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
    M.setZero();
    for(int i=0;i<3;++i)
    {
        M.block(4*i,4*i,4,4)=Eigen::Matrix4d::Ones();
    }
    M*=1.0/10.0*density*volume;
   
 }