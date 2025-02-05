#include "AugmentedLagrangian.hpp"
#include <fstream>

void write_solution(const std::string& file_name, const SolutionType& sol){
    std::ofstream file(file_name);
    file << "Number of iterations: " << sol.final_iter << endl << endl;
    file << "Convergence criterion: " << sol.conv_criterion << endl << endl;
    file << "Result: " << endl << sol.result;
}

int main(){

    MatrixXd A(3,2);
    A << 1, 2, 1, 3, 2, 5;

    MatrixXd B(1,2);
    B << -1, 1;

    VectorXd d(3);
    d << 1, 1, 1;

    VectorXd g(1);
    g << 0;

    double rho = 1;
    double tol = 1e-5;
    int max_iter = 100;

    // Parameters param;
    // param.g0 = 0;
    // param.g1 = 1;
    // param.f = 1;
    // param.rho_hat = 0.5;
    // size_t n = 16;

    MultipliersMethod mm(A, B, g, d, rho, linearSystemSolverType::ColPivHouseholderQR);
    // MultipliersMethod mm(param, n, linearSystemSolverType::ColPivHouseholderQR);

    SolutionType sol = mm.solve(tol, max_iter);

    write_solution("sol.txt", sol);

}