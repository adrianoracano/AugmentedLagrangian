#include "DifferentialProblemSolver.hpp"

//=========================================================================================
// DIFFERENTIAL PROBLEM SOLVER //==========================================================
//=========================================================================================

/* // void DifferentialProblemSolver::setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver)
// {
//     switch (linear_sys_solver) {
//         case linearSystemSolverType::PartialPivLU:
//             linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.partialPivLu());
//             break;
//         case linearSystemSolverType::FullPivLU:
//             linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.fullPivLu());
//             break;
//         case linearSystemSolverType::HouseholderQR:
//             linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.householderQr());
//             break;
//         case linearSystemSolverType::ColPivHouseholderQR:
//             linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.colPivHouseholderQr());
//             break;
//         case linearSystemSolverType::FullPivHouseholderQR:
//             linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.fullPivHouseholderQr());
//             break;
//         case linearSystemSolverType::CompleteOrthogonalDecomposition:
//             linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.completeOrthogonalDecomposition());
//             break;
//         case linearSystemSolverType::LLT:
//             linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.llt());
//             break;
//         case linearSystemSolverType::LDLT:
//             linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.ldlt());
//             break;
//         default:
//             throw std::invalid_argument("Invalid solver type");
//     }
// }
 */

void DifferentialProblemSolver::set_params(const Parameters& params){
        g.resize(2);
        g(0) = params.g0;
        g(1) = params.g1;
        rho_hat = params.rho_hat;
    };

void DifferentialProblemSolver::update_A(){
    A.resize(n, n);
    std::vector<Triplet> nodes;
    nodes.emplace_back(0, 0, 1);
    nodes.emplace_back(0, 1, -1);
    nodes.emplace_back(n - 1, n - 1, 1);
    nodes.emplace_back(n - 1, n - 2, -1);
    for (size_t i = 1; i < n - 1; ++i){
        nodes.emplace_back(i, i-1, -1);
        nodes.emplace_back(i, i, 2);
        nodes.emplace_back(i, i+1, -1);
    }
    A.setFromTriplets(nodes.begin(), nodes.end());
    A.makeCompressed();
};

void DifferentialProblemSolver::update_B(){
    B.resize(2, n);
    std::vector<Triplet> B_coeff;
    B_coeff.emplace_back(0, 0, 1);
    B_coeff.emplace_back(1, n-1, 1);
    B.setFromTriplets(B_coeff.begin(), B_coeff.end());
    B.makeCompressed();
};

void DifferentialProblemSolver::update_f(){
    f.resize(n);
    f.fill(h);
    f(0) /= 2;
    f(n-1) /= 2;
};

void DifferentialProblemSolver::set_n_and_update(const size_t& n_){
    n = n_;
    h = 1.0/(n+1);
    rho = rho_hat / h;
    update_A();
    update_B();
    update_f();
};

DifferentialProblemSolver::DifferentialProblemSolver(const Parameters& params, const double& n_){
    cout << "ok1" << endl;
    set_params(params);
    set_n_and_update(n_);
    LeftMatrix = (A + rho * B.transpose() * B).pruned();
    RightFirstAddendum = f + rho * B.transpose() * g;
    linearSystemSolver = std::make_unique<Eigen::FullPivLU<MatrixXd>>(LeftMatrix.fullPivLu());
    cout << "ok2" << endl;
};

// DifferentialProblemSolver::DifferentialProblemSolver(
//                        const Parameters& params, 
//                        const double& n_,
//                        linearSystemSolverType const& linear_sys_solver): 
//                        DifferentialProblemSolver::DifferentialProblemSolver(params, n_)
// {
//     // setLinearSystemSolver(linear_sys_solver);
// };

// void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver) override;

VectorXd DifferentialProblemSolver::solve(const VectorXd& lambda) const {
    cout << "lhs is  = " << endl << LeftMatrix << endl;
    cout << endl << "rhs is : " << endl << RightFirstAddendum - B.transpose() * lambda << endl;
    // return std::visit([&](auto& solverInstance) -> VectorXd {
    //     cout << endl << "From inside the std::visit lambda : " << endl;
    //     cout << endl << "lhs is  = " << endl << LeftMatrix << endl << endl;
    //     cout << endl << "rhs is : " << endl << RightFirstAddendum - B.transpose() * lambda << endl << endl;
    //     return solverInstance.solve(RightFirstAddendum - B.transpose() * lambda);
    // }, *linearSystemSolver);
    return linearSystemSolver->solve(RightFirstAddendum - B.transpose() * lambda);
}

VectorXd DifferentialProblemSolver::update(const VectorXd& lambda, const VectorXd& x) const {
    return  lambda - rho * (g - B * x);
};

double DifferentialProblemSolver::get_residual_norm(const VectorXd& x) const  {
    return  (A*x - f).norm();
};

double DifferentialProblemSolver::conv_criterion(const VectorXd& x, const VectorXd& x_old) const  {
    return (x - x_old).transpose() * (A + h * MatrixXd::Identity(n, n)) * (x - x_old);
};

