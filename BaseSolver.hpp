/**
 * @file BaseSolver.hpp
 * @brief Base class for all solver implementations
 * @details Provides interface and common functionality for linear system solvers
 */
#ifndef BASESOLVER_HPP
#define BASESOLVER_HPP

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>

#include <iostream>
#include <memory>
#include <string>
#include <functional>
#include <type_traits>
#include <utility>
#include <variant>

using Eigen::MatrixXd;
using SparseMatrixXd = Eigen::SparseMatrix<double>; // compressed sparse row
using Eigen::VectorXd;
using Triplet = Eigen::Triplet<double>;

using std::size_t;
using std::cout; 
using std::cerr; 
using std::endl;

/** @brief Available linear system solver types */
enum class linearSystemSolverType {
    PartialPivLU,              ///< Partial pivoting LU
    FullPivLU,                 ///< Full pivoting LU
    HouseholderQR,            ///< Standard Householder QR
    ColPivHouseholderQR,      ///< Column pivoting QR
    FullPivHouseholderQR,     ///< Full pivoting QR
    CompleteOrthogonalDecomposition,  ///< COD decomposition
    LLT,                      ///< Cholesky decomposition
    LDLT                      ///< LDLT decomposition
};

using SolverVariant = std::variant<
    Eigen::PartialPivLU<MatrixXd>,
    Eigen::FullPivLU<MatrixXd>,
    Eigen::HouseholderQR<MatrixXd>,
    Eigen::ColPivHouseholderQR<MatrixXd>,
    Eigen::FullPivHouseholderQR<MatrixXd>,
    Eigen::CompleteOrthogonalDecomposition<MatrixXd>,
    Eigen::LLT<MatrixXd>,
    Eigen::LDLT<MatrixXd>>;

/**
 * @class SolutionType
 * @brief Container for solver results
 */
class SolutionType {
public:
    /**
     * @brief Constructor for solution container
     * @param final_iter_ Final iteration count
     * @param conv_criterion_ Convergence criterion value
     * @param residual_norm_ Norm of the residual
     * @param result_ Solution vector
     */
    SolutionType(const int& final_iter_, const double& conv_criterion_, const double& residual_norm_, const VectorXd& result_):
    final_iter(final_iter_), conv_criterion(conv_criterion_), result(result_), residual_norm(residual_norm_) {};

    int final_iter;           ///< Number of iterations performed
    double conv_criterion;    ///< Final convergence criterion value
    double residual_norm;     ///< Final residual norm
    VectorXd result;         ///< Solution vector
};

/**
 * @class BaseSolver
 * @brief Abstract base class for all solvers
 * @details Provides interface for solving linear systems with various methods
 */
class BaseSolver {
public:
    BaseSolver() = default;
    virtual ~BaseSolver() {};

    /**
     * @brief Solves the system with given Lagrange multiplier
     * @param lambda Current Lagrange multiplier
     * @return Solution vector
     */
    virtual VectorXd solve(const VectorXd& lambda) const = 0;

    /**
     * @brief Updates solution with new multiplier
     * @param lambda New Lagrange multiplier
     * @param x Current solution
     * @return Updated solution
     */
    virtual VectorXd update(const VectorXd& lambda, const VectorXd& x) const = 0;

    /**
     * @brief Computes residual norm
     * @param x Current solution
     * @return Norm of the residual
     */
    virtual double get_residual_norm(const VectorXd& x) const = 0;

    /**
     * @brief Evaluates convergence criterion
     * @param x Current solution
     * @param x_old Previous solution
     * @return Convergence measure
     */
    virtual double conv_criterion(const VectorXd& x, const VectorXd& x_old) const = 0;

    /**
     * @brief Sets the linear system solver type
     * @param linear_sys_solver Solver type to use
     * @throws std::invalid_argument if solver type is invalid
     */
    void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver);

protected:
    std::unique_ptr<SolverVariant> linearSystemSolver;  ///< Pointer to solver implementation
    MatrixXd LeftMatrix;                               ///< System matrix
};

/**
 * @struct Parameters
 * @brief Problem parameters container
 */
struct Parameters {
    double g0;        ///< First boundary condition
    double g1;        ///< Second boundary condition
    double rho_hat;   ///< Penalty parameter
};

inline void BaseSolver::setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver)
{
    switch (linear_sys_solver) {
        case linearSystemSolverType::PartialPivLU:
            linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.partialPivLu());
            break;
        case linearSystemSolverType::FullPivLU:
            linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.fullPivLu());
            break;
        case linearSystemSolverType::HouseholderQR:
            linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.householderQr());
            break;
        case linearSystemSolverType::ColPivHouseholderQR:
            linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.colPivHouseholderQr());
            break;
        case linearSystemSolverType::FullPivHouseholderQR:
            linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.fullPivHouseholderQr());
            break;
        case linearSystemSolverType::CompleteOrthogonalDecomposition:
            linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.completeOrthogonalDecomposition());
            break;
        case linearSystemSolverType::LLT:
            linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.llt());
            break;
        case linearSystemSolverType::LDLT:
            linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.ldlt());
            break;
        default:
            throw std::invalid_argument("Invalid solver type");
    }
}

#endif // BASESOLVER_HPP