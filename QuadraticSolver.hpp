#ifndef QUADRATICSOLVER_HPP
#define QUADRATICSOLVER_HPP

#include "BaseSolver.hpp"

class QuadraticSolver : public BaseSolver {

    public:

    QuadraticSolver() {};
    QuadraticSolver(MatrixXd A_, MatrixXd B_, VectorXd g_, VectorXd d_, double rho_);
    QuadraticSolver(MatrixXd A_, MatrixXd B_, VectorXd g_, VectorXd d_, double rho_, linearSystemSolverType const& linear_sys_solver);

    // void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver);
    
    void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver);
    virtual VectorXd solve(const VectorXd& lambda) const override;
    virtual VectorXd update(const VectorXd& lambda, const VectorXd& x) const override;
    virtual double get_residual_norm(const VectorXd& x) const override;
    virtual double conv_criterion(const VectorXd& x, const VectorXd& x_old) const override;
    
    virtual ~QuadraticSolver() {};

    private:
    MatrixXd LeftMatrix;
    VectorXd RightFirstAddendum;
    std::unique_ptr<SolverVariant> linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.partialPivLu());
    MatrixXd A, B;
    VectorXd g, d;
    double rho;
};

#endif //QUADRATICSOLVER_HPP

