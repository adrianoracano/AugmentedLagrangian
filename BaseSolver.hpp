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


enum class linearSystemSolverType {PartialPivLU, 
                          FullPivLU,
                          HouseholderQR,
                          ColPivHouseholderQR,
                          FullPivHouseholderQR,
                          CompleteOrthogonalDecomposition,
                          LLT,
                          LDLT};

using SolverVariant = std::variant<
    Eigen::PartialPivLU<MatrixXd>,
    Eigen::FullPivLU<MatrixXd>,
    Eigen::HouseholderQR<MatrixXd>,
    Eigen::ColPivHouseholderQR<MatrixXd>,
    Eigen::FullPivHouseholderQR<MatrixXd>,
    Eigen::CompleteOrthogonalDecomposition<MatrixXd>,
    Eigen::LLT<MatrixXd>,
    Eigen::LDLT<MatrixXd>>;

class SolutionType{
    public:
    SolutionType(const int& final_iter_, const double& conv_criterion_, const double& residual_norm_, const VectorXd& result_):
    final_iter(final_iter_), conv_criterion(conv_criterion_), result(result_), residual_norm(residual_norm_) {};

    int final_iter;
    double conv_criterion;
    double residual_norm;
    VectorXd result;
};

class BaseSolver{
    public:
    BaseSolver() = default;
    virtual ~BaseSolver() {};
    virtual VectorXd solve(const VectorXd& lambda) const = 0;
    virtual VectorXd update(const VectorXd& lambda, const VectorXd& x) const = 0;
    virtual double get_residual_norm(const VectorXd& x) const = 0;
    virtual double conv_criterion(const VectorXd& x, const VectorXd& x_old) const = 0;
    // virtual void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver);

    // std::unique_ptr<SolverVariant> linearSystemSolver;
    // MatrixXd LeftMatrix;
};

struct Parameters{
    double g0;
    double g1;
    double rho_hat;
};

#endif // BASESOLVER_HPP