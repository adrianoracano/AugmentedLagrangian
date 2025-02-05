#ifndef METHODOFMULTIPLIERS_HPP
#define METHODOFMULTIPLIERS_HPP

#include "QuadraticSolver.hpp"
#include "DifferentialProblemSolver.hpp"

class MultipliersMethod{

    public:  
    MultipliersMethod(MatrixXd A, MatrixXd B, VectorXd g, VectorXd d_vec, double rho, 
                      linearSystemSolverType const& linear_sys_solver = linearSystemSolverType::FullPivLU): m(A.rows()), n(A.cols()), d(B.rows()){
        
        solver = std::make_unique<QuadraticSolver>(A, B, g, d_vec, rho, linear_sys_solver);
    };

    MultipliersMethod(const Parameters& params, const double& n_,
                      linearSystemSolverType const& linear_sys_solver = linearSystemSolverType::FullPivLU): m(n_), n(n_), d(2){
        // solver = std::make_unique<DifferentialProblemSolver>(params, n);
        cout << "Entered constructor of MultipliersMethod" << endl;
        solver = std::make_unique<DifferentialProblemSolver>(params, n);
    };

    SolutionType solve(double tol, int max_iter);
    
    private:
    std::unique_ptr<BaseSolver> solver;
    // std::unique_ptr<DifferentialProblemSolver> solver;
    int m, n, d;
};

#endif //METHODOFMULTIPLIERS_HPP
