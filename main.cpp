// #include "AugmentedLagrangian.hpp"
#include "MethodOfMultipliers.hpp"
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

    MultipliersMethod mm(A, B, g, d, rho, linearSystemSolverType::ColPivHouseholderQR);

    SolutionType sol = mm.solve(tol, max_iter);

    write_solution("sol.txt", sol);

}