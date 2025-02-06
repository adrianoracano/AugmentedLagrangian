/**
 * @file MethodOfMultipliers.hpp
 * @brief Implementation of the Method of Multipliers (Augmented Lagrangian method)
 * @details Provides solution methods for both quadratic and differential problems using
 *          the method of multipliers optimization approach.
 */

#ifndef METHODOFMULTIPLIERS_HPP
#define METHODOFMULTIPLIERS_HPP

#include "QuadraticSolver.hpp"
#include "DifferentialProblemSolver.hpp"

/**
 * @class MultipliersMethod
 * @brief Implements the method of multipliers optimization algorithm
 * 
 * @details This class provides a unified interface for solving constrained optimization
 *          problems using either QuadraticSolver or DifferentialProblemSolver.
 */
class MultipliersMethod {
public:  
    /**
     * @brief Constructor for quadratic programming problems
     * @param A System matrix for quadratic term
     * @param B Constraint matrix
     * @param g Linear term vector
     * @param d_vec Right-hand side vector for constraints
     * @param rho Penalty parameter
     * @param linear_sys_solver Type of linear system solver to use
     */
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
    std::unique_ptr<BaseSolver> solver;  ///< Pointer to solver implementation
    int m;                               ///< Number of constraints
    int n;                               ///< Problem dimension
    int d;                               ///< Degrees of freedom
};

#endif //METHODOFMULTIPLIERS_HPP
