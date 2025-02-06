// #endif //DIFFERENTIALPROBLEMSOLVER_HPP

/**
 * @file DifferentialProblemSolver.hpp
 * @brief Solver for differential problems in augmented Lagrangian framework
 */

#ifndef DIFFERENTIALPROBLEMSOLVER_HPP
#define DIFFERENTIALPROBLEMSOLVER_HPP

#include "BaseSolver.hpp"

/**
 * @class DifferentialProblemSolver
 * @brief Implements solver for differential equations using augmented Lagrangian method
 * 
 * Inherits from BaseSolver to provide specific implementation for differential problems.
 * Uses sparse matrices for efficient computation of large systems.
 */
class DifferentialProblemSolver : public BaseSolver {
public:
    /** @brief Default constructor */
    DifferentialProblemSolver() {};

    /**
     * @brief Constructor with parameters
     * @param params Problem parameters
     * @param n_ Number of discretization points
     */
    DifferentialProblemSolver(const Parameters& params, const double& n_);

    /**
     * @brief Full constructor with linear system solver
     * @param params Problem parameters
     * @param n_ Number of discretization points
     * @param linear_sys_solver Linear system solver implementation
     */
    DifferentialProblemSolver(const Parameters& params, 
                             const double& n_,
                             linearSystemSolverType const& linear_sys_solver);

    /**
     * @brief Updates problem parameters
     * @param params New parameters
     */
    void set_params(const Parameters& params);

    /** @brief Updates matrix A in the system */
    void update_A();

    /** @brief Updates matrix B in the system */
    void update_B();

    /** @brief Updates right-hand side vector f */
    void update_f();

    /**
     * @brief Updates discretization and rebuilds system
     * @param n_ New number of discretization points
     */
    void set_n_and_update(const size_t& n_);

    /**
     * @brief Solves the differential problem
     * @param lambda Lagrange multiplier vector
     * @return Solution vector
     */
    VectorXd solve(const VectorXd& lambda) const override;

    /**
     * @brief Updates solution with new multiplier
     * @param lambda New Lagrange multiplier
     * @param x Current solution
     * @return Updated solution
     */
    VectorXd update(const VectorXd& lambda, const VectorXd& x) const override;

    /**
     * @brief Computes residual norm
     * @param x Current solution
     * @return Norm of the residual
     */
    double get_residual_norm(const VectorXd& x) const override;

    /**
     * @brief Evaluates convergence criterion
     * @param x Current solution
     * @param x_old Previous solution
     * @return Convergence measure
     */
    double conv_criterion(const VectorXd& x, const VectorXd& x_old) const override;

private:
    using BaseSolver::setLinearSystemSolver;
    using BaseSolver::linearSystemSolver;
    using BaseSolver::LeftMatrix;

    VectorXd RightFirstAddendum;     ///< First term of right-hand side
    SparseMatrixXd A, B;             ///< System matrices
    size_t n;                        ///< Number of discretization points
    VectorXd g, f;                   ///< System vectors
    double rho_hat;                  ///< Penalty parameter
    double rho;                      ///< Augmented Lagrangian parameter
    double h;                        ///< Grid spacing
};

#endif //DIFFERENTIALPROBLEMSOLVER_HPP