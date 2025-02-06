/**
 * @file QuadraticSolver.hpp
 * @brief Solver for quadratic optimization problems
 * @details Implements solution methods for problems of the form:
 *          min 1/2 x^T A x + g^T x
 *          s.t. Bx = d
 */
#ifndef QUADRATICSOLVER_HPP
#define QUADRATICSOLVER_HPP

#include "BaseSolver.hpp"

/**
 * @class QuadraticSolver
 * @brief Implements solver for quadratic optimization problems
 * @details Inherits from BaseSolver to provide specific implementation
 *          for quadratic programming problems
 */
class QuadraticSolver : public BaseSolver {
public:
    /** @brief Default constructor */
    QuadraticSolver() {};

    /**
     * @brief Constructor with problem data
     * @param A_ System matrix for quadratic term
     * @param B_ Constraint matrix
     * @param g_ Linear term vector
     * @param d_ Right-hand side vector
     * @param rho_ Penalty parameter
     */
    QuadraticSolver(MatrixXd A_, MatrixXd B_, VectorXd g_, VectorXd d_, double rho_);

    /**
     * @brief Full constructor with solver specification
     * @param A_ System matrix for quadratic term
     * @param B_ Constraint matrix
     * @param g_ Linear term vector
     * @param d_ Right-hand side vector
     * @param rho_ Penalty parameter
     * @param linear_sys_solver Linear system solver type
     */
    QuadraticSolver(MatrixXd A_, MatrixXd B_, VectorXd g_, VectorXd d_, double rho_, 
                    linearSystemSolverType const& linear_sys_solver);

    virtual VectorXd solve(const VectorXd& lambda) const override;
    virtual VectorXd update(const VectorXd& lambda, const VectorXd& x) const override;
    virtual double get_residual_norm(const VectorXd& x) const override;
    virtual double conv_criterion(const VectorXd& x, const VectorXd& x_old) const override;
    
    virtual ~QuadraticSolver() {};

private:
    using BaseSolver::setLinearSystemSolver;   ///< Use base class solver setter
    using BaseSolver::linearSystemSolver;       ///< Use base class solver
    using BaseSolver::LeftMatrix;
    VectorXd RightFirstAddendum;
    MatrixXd A, B;
    VectorXd g, d;
    double rho;
};

#endif //QUADRATICSOLVER_HPP

