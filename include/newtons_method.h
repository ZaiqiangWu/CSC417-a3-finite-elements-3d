#include <Eigen/Dense>
#include <EigenTypes.h>
#include "timer.h"
#define show_time

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {
   Eigen::VectorXd d;
    Eigen::VectorXd new_tmp_g;
    Eigen::VectorXd tmp_x;

   double c1=1e-4;
   double c2=0.9;
   //x0 is velocity!!!
   Timer timer=Timer();


    for(int i=0;i<maxSteps;++i)
   {
        timer.begin();
        g(tmp_g,x0);
#ifdef show_time
        std::cout<<"time for evaluating gradient: "<<timer.get_elapsed_time()<<std::endl;
#endif

        timer.begin();
        H(tmp_H,x0);
#ifdef show_time
        std::cout<<"time for evaluating Hessian: "<<timer.get_elapsed_time()<<std::endl;
#endif

        timer.begin();
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(tmp_H);
#ifdef show_time
        std::cout<<"time for factorize H (only for direct solver): "<<timer.get_elapsed_time()<<std::endl;
#endif
        timer.begin();
        d=solver.solve(-tmp_g);
        //d=-tmp_g;

#ifdef show_time
        std::cout<<"time for back substitution (only for direct solver): "<<timer.get_elapsed_time()<<std::endl;
#endif
        //line search
        double alpha=1.0;
        for(int j=0;j<10;++j)
        {
            tmp_x=x0+alpha*d;
            g(new_tmp_g,tmp_x);
            if(f(tmp_x)<f(x0)+c1*alpha*d.dot(tmp_g)&&(-d.dot(new_tmp_g)<-c2*d.dot(tmp_g)))
            //if(f(tmp_x)<f(x0))
                break;
            else
            {
                alpha=alpha*0.5;
                std::cout<<"step too large:"<<j<<std::endl;
            }

        }
        x0=tmp_x;
   }




   return 0.0;
}
