#include "DifferentialProblemSolver.hpp"

//=========================================================================================
// DIFFERENTIAL PROBLEM SOLVER //==========================================================
//=========================================================================================

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
    A /= h;
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
    h = 1.0/(n-1);
    rho = rho_hat / h;
    update_A();
    update_B();
    update_f();
};

DifferentialProblemSolver::DifferentialProblemSolver(const Parameters& params, const double& n_){
    // cout << "ok1" << endl;
    set_params(params);
    set_n_and_update(n_);
    LeftMatrix = (A + rho * B.transpose() * B).pruned();
    RightFirstAddendum = f + rho * B.transpose() * g;
    linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.fullPivLu());
    // cout << "ok2" << endl;
};

DifferentialProblemSolver::DifferentialProblemSolver(
                       const Parameters& params, 
                       const double& n_,
                       linearSystemSolverType const& linear_sys_solver): 
                       DifferentialProblemSolver::DifferentialProblemSolver(params, n_)
{
    setLinearSystemSolver(linear_sys_solver);
};

// void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver) override;

VectorXd DifferentialProblemSolver::solve(const VectorXd& lambda) const {
    return std::visit([&](auto& solverInstance) -> VectorXd {
        return solverInstance.solve(RightFirstAddendum - B.transpose() * lambda);
    }, *linearSystemSolver);
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

