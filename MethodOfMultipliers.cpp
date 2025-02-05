#include "MethodOfMultipliers.hpp"

//=========================================================================================
// METHOD OF MULTIPLIERS //================================================================
//=========================================================================================

SolutionType MultipliersMethod::solve(double tol, int max_iter){
    cout << "using MultipliersMethod::solve" << endl;
    VectorXd lambda = VectorXd::Zero(d);
    VectorXd lambda_new;
    VectorXd x_old, x{VectorXd::Zero(n)};
    size_t it = 0;
    double err = 2 * tol;

    while (it < max_iter && err > tol){
        x_old = x;
        cout << "iteration n " << it << endl;
        x = solver->solve(lambda);
        lambda_new = solver->update(lambda, x);
        err = solver->conv_criterion(x, x_old);
        lambda = lambda_new;
        it++;
    }
    if (err >= tol)
        cerr << "Not converged, final iterations : " << it << endl;
    cout << "Final iterations : " << it << endl;
    cout << "Solution found : " << endl << x << endl;
    return SolutionType(it, err, solver->get_residual_norm(x), x);
}








