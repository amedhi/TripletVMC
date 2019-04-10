#include <Eigen/Core>
#include <iostream>
#include <iomanip>
#include "./cppoptlib/meta.h"
#include "./cppoptlib/boundedproblem.h"
//#include "./cppoptlib/solver/bfgssolver.h"
//#include "./cppoptlib/solver/conjugatedgradientdescentsolver.h"
//#include "./cppoptlib/solver/newtondescentsolver.h"
//#include "./cppoptlib/solver/neldermeadsolver.h"
//#include "./cppoptlib/solver/lbfgssolver.h"
#include "./cppoptlib/solver/lbfgsbsolver.h"
//#include "./cppoptlib/solver/cmaessolver.h"


using Eigen::VectorXf;
using Eigen::MatrixXf;
//using namespace LBFGSpp;
using namespace cppoptlib; 

class Rosenbrock : public BoundedProblem<double> {
  public:
    using TVector = Problem<double>::TVector;
    double operator()(const TVector &x, TVector &grad)
    {
       // gradient
        grad[0]  = -2 * (1 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2 * x[0]);
        grad[1]  =                   200 * (x[1] - x[0] * x[0]);

        // function value
        const double t1 = (1 - x[0]);
        const double t2 = (x[1] - x[0] * x[0]);
        return   t1 * t1 + 100 * t2 * t2;
    }

    double operator()(const TVector &x)
    {
        return value(x);
    }
    // this is just the objective (NOT optional)
    double value(const TVector &x) {
        const double t1 = (1 - x[0]);
        const double t2 = (x[1] - x[0] * x[0]);
        return   t1 * t1 + 100 * t2 * t2;
    }

    bool callback(const cppoptlib::Criteria<double> &state, const TVector &x) {
        std::cout << "(" << std::setw(2) << state.iterations << ")"
                  << " ||dx|| = " << std::fixed << std::setw(8) << std::setprecision(4) << state.gradNorm
                  << " ||x|| = "  << std::setw(6) << x.norm()
                  << " f(x) = "   << std::setw(8) << value(x)
                  << " x = [" << std::setprecision(8) << x.transpose() << "]" << std::endl;
        return true;
    }
};

/*
class Rosenbrock
{
private:
    int n;
public:
    Rosenbrock(int n_) : n(n_) {}
    float operator()(const VectorXf& x, VectorXf& grad)
    {
        float fx = 0.0;
        for(int i = 0; i < n; i += 2)
        {
            float t1 = 1.0 - x[i];
            float t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
            fx += t1 * t1 + t2 * t2;
        }
        return fx;
    }
};
*/


int main()
{
    /*
    const int n = 10;
    LBFGSParam<float> param;
    LBFGSSolver<float> solver(param);
    Rosenbrock fun(n);
    VectorXf x = VectorXf::Zero(n);
    float fx;
    int niter = solver.minimize(fun, x, fx);
    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    */

    Rosenbrock func;
    //func.setDimension(2);
    // choose a starting point
    Eigen::VectorXd x(2); x << -1, 2;
    Eigen::VectorXd lower(2); lower << -1, -1;
    Eigen::VectorXd upper(2); upper << +2, +2;
    func.setBoxConstraint(lower,upper);
  
    LbfgsbSolver<Rosenbrock> solver;
    solver.minimize(func, x);

    std::cout << "argmin      " << x.transpose() << std::endl;
    std::cout << "f in argmin " << func(x) << std::endl;

    return 0;
}
