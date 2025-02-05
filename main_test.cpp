#include "MethodOfMultipliers.hpp"
#include <fstream>

#include <functional>
#include <algorithm>
#include <gnuplot-iostream.hpp>

double u_ex(const double& x){
    return 1/2 * (3*x - x*x);
}

double L2_error(const VectorXd& sol, std::function<double(const double&)> u_ex, size_t N){

    std::vector<double> u(sol.data(), sol.data() + sol.size());
    
    // Piecewise Linear interpolant function of u using std::lerp
    auto linear_interpolant = [&u](double x) -> double {
        // Ensure x is within the bounds [0, 1]
        if (x < 0.0 || x > 1.0) {
            throw std::out_of_range("x is out of the interpolation range [0, 1]");
        }

        // Find the interval [i, i+1]
        int n = u.size() - 1;
        double scaled_x = x * n; // Scale x to match the size of u
        int i = static_cast<int>(std::floor(scaled_x));
        int i1 = std::min(i + 1, n);

        // Compute the interpolation factor t
        double t = scaled_x - i;

        // Perform linear interpolation using std::lerp
        double y = std::lerp(u[i], u[i1], t);
        return y;
    };
    
    double error = 0;

    for (size_t i = 0; i < N; ++i){
        double x = static_cast<double>(i) / N;
        double u_ex_val = u_ex(x);
        double u_approx_val = linear_interpolant(x);
        error += (u_ex_val - u_approx_val) * (u_ex_val - u_approx_val);
    }

    return sqrt(1.0 / N * error);
}

void write_solution(const std::string& file_name, const SolutionType& sol, size_t N){
    std::ofstream file(file_name);
    file << "Number of iterations: " << sol.final_iter << endl << endl;
    file << "Convergence criterion: " << sol.conv_criterion << endl << endl;
    file << "Result: " << endl << sol.result << endl << endl;
    file << "Error: " << L2_error(sol.result, u_ex, N);
}

int main(){

    double tol = 1e-5;
    int max_iter = 100;

    Parameters param;
    param.g0 = 0;
    param.g1 = 1;
    param.rho_hat = 0.5;
    size_t n, n_evaluation_points;

    MultipliersMethod mm(param, 16, linearSystemSolverType::ColPivHouseholderQR);
    SolutionType sol = mm.solve(tol, max_iter);
    // write_solution("prova.txt", sol, 16 * 10);

    // Gnuplot gp; // Create a Gnuplot object

    // // Specify terminal type for better window management
    // gp << "set terminal wxt size 800,600 enhanced\n";

    // // Prepare a vector of plot data and titles for Gnuplot
    // std::vector<std::vector<std::pair<double, double>>> all_plot_data;
    // std::vector<std::string> titles;

    // for (size_t i = 4; i <= 6; i++) {
    //     n = std::pow(2, i);
    //     cout << "n = " << n << endl;
    //     MultipliersMethod mm(param, n, linearSystemSolverType::ColPivHouseholderQR);
    //     SolutionType sol = mm.solve(tol, max_iter);

    //     n_evaluation_points = n * 10;
    //     std::string filename = "sol_n_" + std::to_string(n);
    //     write_solution(filename, sol, n_evaluation_points);

    //     // Prepare data for plotting
    //     std::vector<std::pair<double, double>> plot_data;
    //     for (size_t j = 0; j < n; ++j) {
    //         double x = static_cast<double>(j) / n;
    //         plot_data.push_back({x, sol.result[j]});
    //     }

    //     all_plot_data.push_back(plot_data);
    //     titles.push_back("n = " + std::to_string(n));
    // }

    // // Plot all solutions in the same window with legends
    // gp << "set title 'Solutions for different n'\n";
    // gp << "plot";
    // for (size_t i = 0; i < all_plot_data.size(); ++i) {
    //     if (i > 0) gp << ", ";
    //     gp << "'-' with lines title '" << titles[i] << "'";
    // }
    // gp << "\n";
    // for (const auto& plot_data : all_plot_data) {
    //     gp.send1d(plot_data);
    // }

    // // Flush the commands to ensure they are all processed
    // gp.flush();

    // // Keep the window open until the user closes it
    // std::cout << "Press enter to close the plot window...";
    // std::cin.get();

    return 0;
}