#include <json.hpp>
#include "MethodOfMultipliers.hpp"
#include <fstream>
#include <functional>
#include <algorithm>
#include <gnuplot-iostream.hpp>

using json = nlohmann::json;

double u_ex(const double& x) {
    return 1.0/2 * (3*x - x*x);
}

double L2_error(const VectorXd& sol, std::function<double(const double&)> u_ex, size_t N) {
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

    for (size_t i = 0; i < N; ++i) {
        double x = static_cast<double>(i) / N;
        double u_ex_val = u_ex(x);
        double u_approx_val = linear_interpolant(x);
        error += (u_ex_val - u_approx_val) * (u_ex_val - u_approx_val);
    }

    return sqrt(1.0 / N * error);
}

void write_solution(const std::string& file_name, const SolutionType& sol, size_t N) {
    std::ofstream file(file_name);
    file << "Number of iterations: " << sol.final_iter << std::endl << std::endl;
    file << "Convergence criterion: " << sol.conv_criterion << std::endl << std::endl;
    file << "Result: " << std::endl << sol.result << std::endl << std::endl;
    file << "Error: " << L2_error(sol.result, u_ex, N);
}

int main() {
    // Read JSON configuration
    std::ifstream config_file("data.json");
    json config;
    config_file >> config;

    // Extract parameters
    const double tol = config["solver_params"]["tol"];
    const int max_iter = config["solver_params"]["max_iter"];

    Parameters params;
    // Default values if not specified in data.json
    params.rho_hat = config["solver_params"].value("rho_hat", 0.5);    
    params.g0 = config["solver_params"].value("g0", 0.0);
    params.g1 = config["solver_params"].value("g1", 1.0);              

    const std::vector<int> n_values = config["discretization"]["n_values"];
    const size_t n_exact_points = config["discretization"]["n_exact_points"];
    const std::string base_filename = config["output"]["base_filename"];
    const std::string extension = config["output"]["extension"];

    // Initialize plotting data
    std::vector<std::vector<std::pair<double, double>>> all_plot_data;
    std::vector<std::string> titles;

    // Solve for different n values
    for (size_t n : n_values) {

        MultipliersMethod mm(params, n);
        auto sol = mm.solve(tol, max_iter);

        // Write solution to file
        std::string filename = base_filename + std::to_string(n) + extension;
        write_solution(filename, sol, n);

        // Store plot data
        std::vector<std::pair<double, double>> plot_data;
        for (size_t j = 0; j < n; ++j) {
            double x = static_cast<double>(j * (n)) / (n * (n-1));
            plot_data.push_back({x, sol.result[j]});
        }

        all_plot_data.push_back(plot_data);
        titles.push_back("n = " + std::to_string(n));
    }

    // Prepare data for plotting
    std::vector<std::pair<double, double>> exact_sol_data;
    for (size_t j = 0; j <= n_exact_points; ++j) {
        double x = static_cast<double>(j) / (n_exact_points);
        exact_sol_data.push_back({x, u_ex(x)});
    }

    // Plot all solutions in the same window with legends
    Gnuplot gp; // Create a Gnuplot object

    // Specify terminal type for better window management
    gp << "set terminal wxt size 800,600 enhanced\n";
    gp << "set title '" << config["plot"]["title"] << "'\n";
    gp << "plot";
    for (size_t i = 0; i < all_plot_data.size(); ++i) {
        if (i > 0) gp << ", ";
        gp << "'-' with lines title '" << titles[i] << "'";
    }

    // Add the exact solution plot
    gp << ", '-' with linespoints title 'Exact Solution' linewidth 2 pointtype 7 pointsize 1.5\n";

    // Send data to Gnuplot
    for (const auto& plot_data : all_plot_data) {
        gp.send1d(plot_data);
    }
    gp.send1d(exact_sol_data);

    // Flush the commands to ensure they are all processed
    gp.flush();

    // Keep the window open until the user closes it
    std::cout << "Press enter to close the plot window...";
    std::cin.get();

    return 0;
}