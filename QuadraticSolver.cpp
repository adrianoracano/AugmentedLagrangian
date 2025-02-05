#include "QuadraticSolver.hpp"

//=========================================================================================
// QUADRATIC SOLVER //=====================================================================
//=========================================================================================

void QuadraticSolver::setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver)
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
            cout << "Using this solver" << endl;
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


QuadraticSolver::QuadraticSolver(MatrixXd A_, MatrixXd B_, VectorXd g_, VectorXd d_, double rho_) :rho(rho_), A(A_), B(B_), d(d_), g(g_)
    {
        LeftMatrix = A.transpose() * A + rho * B.transpose() * B;
        RightFirstAddendum = A.transpose() * d + rho * B.transpose() * g;
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