#ifndef DIFFERENTIALPROBLEMSOLVER_HPP
#define DIFFERENTIALPROBLEMSOLVER_HPP

#include "BaseSolver.hpp"

class DifferentialProblemSolver : public BaseSolver {

    public:
    DifferentialProblemSolver() {};
    DifferentialProblemSolver(const Parameters& params, 
                       const double& n_);
    // DifferentialProblemSolver(
    //                    const Parameters& params, 
    //                    const double& n_,
    //                    linearSystemSolverType const& linear_sys_solver);

    void set_params(const Parameters& params);
    void update_A();
    void update_B();
    void update_f();
    void set_n_and_update(const size_t& n_);

    // void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver) override;
    
    // void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver);

    VectorXd solve(const VectorXd& lambda) const override;
    VectorXd update(const VectorXd& lambda, const VectorXd& x) const override;
    double get_residual_norm(const VectorXd& x) const override;
    double conv_criterion(const VectorXd& x, const VectorXd& x_old) const override;

    private:
    MatrixXd LeftMatrix;
    VectorXd RightFirstAddendum;
    // std::unique_ptr<SolverVariant> linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.partialPivLu());
    std::unique_ptr<Eigen::FullPivLU<MatrixXd>> linearSystemSolver; // = std::make_unique<Eigen::FullPivLU<MatrixXd>>(LeftMatrix.fullPivLu());
    SparseMatrixXd A, B;
    size_t n;
    VectorXd g, f;
    double rho_hat, rho, h;
};


#endif //DIFFERENTIALPROBLEMSOLVER_HPP