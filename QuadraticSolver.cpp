#include "QuadraticSolver.hpp"

//=========================================================================================
// QUADRATIC SOLVER //=====================================================================
//=========================================================================================

QuadraticSolver::QuadraticSolver(MatrixXd A_, MatrixXd B_, VectorXd g_, VectorXd d_, double rho_) :rho(rho_), A(A_), B(B_), d(d_), g(g_)
    {
        LeftMatrix = A.transpose() * A + rho * B.transpose() * B;
        RightFirstAddendum = A.transpose() * d + rho * B.transpose() * g;
        setLinearSystemSolver(linearSystemSolverType::ColPivHouseholderQR);
    }

QuadraticSolver::QuadraticSolver(MatrixXd A_, MatrixXd B_, VectorXd g_, VectorXd d_, double rho_, linearSystemSolverType const& linear_sys_solver):
    QuadraticSolver::QuadraticSolver(A_, B_ , g_, d_, rho_){
        setLinearSystemSolver(linear_sys_solver);
    }
VectorXd QuadraticSolver::solve(const VectorXd& lambda) const {
    return LeftMatrix.colPivHouseholderQr().solve(RightFirstAddendum - B.transpose() *lambda);
};

VectorXd QuadraticSolver::update(const VectorXd& lambda, const VectorXd& x) const {
    return  lambda - rho * (g - B * x);
};

double QuadraticSolver::get_residual_norm(const VectorXd& x) const  {
    return  (A*x - d).norm();
};

double QuadraticSolver::conv_criterion(const VectorXd& x, const VectorXd& x_old) const  {
    return (x - x_old).norm();
};