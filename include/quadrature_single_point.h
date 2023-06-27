#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template<typename Ret, typename Integrand_Func>
inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q, 
                                               Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                               Integrand_Func integrand) {
    Eigen::Vector3d v[4];
    Eigen::Vector3d X=Eigen::Vector3d::Zero();
    for(int i=0;i<4;++i)
    {
        v[i]=q.block(3*element(i),0,3,1);
        X+=v[i];
    }
    X/=4.0;

    integrand(integrated,q,element,X);
    integrated*=volume;

}

