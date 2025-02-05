#ifndef AUGMENTEDLAGRANGIAN_HPP
#define AUGMENTEDLAGRANGIAN_HPP

// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wunused-parameter"

// #include <Eigen/Core>
// #include <Eigen/Dense>
// #include <Eigen/Sparse>
// #include <Eigen/SparseCore>

// #include <iostream>
// #include <memory>
// #include <string>
// #include <functional>
// #include <type_traits>
// #include <utility>
// #include <variant>

// using Eigen::MatrixXd;
// using SparseMatrixXd = Eigen::SparseMatrix<double>; // compressed sparse row
// using Eigen::VectorXd;
// using Triplet = Eigen::Triplet<double>;

// using std::size_t;
// using std::cout; 
// using std::cerr; 
// using std::endl;


// enum class linearSystemSolverType {PartialPivLU, 
//                           FullPivLU,
//                           HouseholderQR,
//                           ColPivHouseholderQR,
//                           FullPivHouseholderQR,
//                           CompleteOrthogonalDecomposition,
//                           LLT,
//                           LDLT};

// using SolverVariant = std::variant<
//     Eigen::PartialPivLU<MatrixXd>,
//     Eigen::FullPivLU<MatrixXd>,
//     Eigen::HouseholderQR<MatrixXd>,
//     Eigen::ColPivHouseholderQR<MatrixXd>,
//     Eigen::FullPivHouseholderQR<MatrixXd>,
//     Eigen::CompleteOrthogonalDecomposition<MatrixXd>,
//     Eigen::LLT<MatrixXd>,
//     Eigen::LDLT<MatrixXd>>;

// class SolutionType{
//     public:
//     SolutionType(const int& final_iter_, const double& conv_criterion_, const double& residual_norm_, const VectorXd& result_):
//     final_iter(final_iter_), conv_criterion(conv_criterion_), result(result_), residual_norm(residual_norm_) {};

//     int final_iter;
//     double conv_criterion;
//     double residual_norm;
//     VectorXd result;
// };

// class BaseSolver{
//     public:
//     BaseSolver() = default;
//     virtual ~BaseSolver() {};
//     virtual VectorXd solve(const VectorXd& lambda) const = 0;
//     virtual VectorXd update(const VectorXd& lambda, const VectorXd& x) const = 0;
//     virtual double get_residual_norm(const VectorXd& x) const = 0;
//     virtual double conv_criterion(const VectorXd& x, const VectorXd& x_old) const = 0;
//     // virtual void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver);

//     // std::unique_ptr<SolverVariant> linearSystemSolver;
//     // MatrixXd LeftMatrix;
// };

// // class QuadraticSolver : public BaseSolver {

// //     public:

// //     QuadraticSolver() {};
// //     QuadraticSolver(MatrixXd A_, MatrixXd B_, VectorXd g_, VectorXd d_, double rho_);
// //     QuadraticSolver(MatrixXd A_, MatrixXd B_, VectorXd g_, VectorXd d_, double rho_, linearSystemSolverType const& linear_sys_solver);

// //     // void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver);
    
// //     void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver);
// //     virtual VectorXd solve(const VectorXd& lambda) const override;
// //     virtual VectorXd update(const VectorXd& lambda, const VectorXd& x) const override;
// //     virtual double get_residual_norm(const VectorXd& x) const override;
// //     virtual double conv_criterion(const VectorXd& x, const VectorXd& x_old) const override;
    
// //     virtual ~QuadraticSolver() {};

// //     private:
// //     MatrixXd LeftMatrix;
// //     VectorXd RightFirstAddendum;
// //     std::unique_ptr<SolverVariant> linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.partialPivLu());
// //     MatrixXd A, B;
// //     VectorXd g, d;
// //     double rho;
// // };

// struct Parameters{
//     double g0;
//     double g1;
//     double rho_hat;
// };

// class DifferentialProblemSolver : public BaseSolver {

//     public:
//     DifferentialProblemSolver() {};
//     DifferentialProblemSolver(const Parameters& params, 
//                        const double& n_);
//     // DifferentialProblemSolver(
//     //                    const Parameters& params, 
//     //                    const double& n_,
//     //                    linearSystemSolverType const& linear_sys_solver);

//     void set_params(const Parameters& params);
//     void update_A();
//     void update_B();
//     void update_f();
//     void set_n_and_update(const size_t& n_);

//     // void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver) override;
    
//     // void setLinearSystemSolver(linearSystemSolverType const& linear_sys_solver);

//     VectorXd solve(const VectorXd& lambda) const override;
//     VectorXd update(const VectorXd& lambda, const VectorXd& x) const override;
//     double get_residual_norm(const VectorXd& x) const override;
//     double conv_criterion(const VectorXd& x, const VectorXd& x_old) const override;

//     private:
//     MatrixXd LeftMatrix;
//     VectorXd RightFirstAddendum;
//     // std::unique_ptr<SolverVariant> linearSystemSolver = std::make_unique<SolverVariant>(LeftMatrix.partialPivLu());
//     std::unique_ptr<Eigen::FullPivLU<MatrixXd>> linearSystemSolver = std::make_unique<Eigen::FullPivLU<MatrixXd>>(LeftMatrix.fullPivLu());
//     SparseMatrixXd A, B;
//     size_t n;
//     VectorXd g, f;
//     double rho_hat, rho, h;
// };

// class MultipliersMethod{

//     public:  
//     MultipliersMethod(MatrixXd A, MatrixXd B, VectorXd g, VectorXd d_vec, double rho, 
//                       linearSystemSolverType const& linear_sys_solver = linearSystemSolverType::FullPivLU): m(A.rows()), n(A.cols()), d(B.rows()){
        
//         solver = std::make_unique<QuadraticSolver>(A, B, g, d_vec, rho, linear_sys_solver);
//     };

//     MultipliersMethod(const Parameters& params, const double& n_,
//                       linearSystemSolverType const& linear_sys_solver = linearSystemSolverType::FullPivLU): m(n_), n(n_), d(2){
//         // solver = std::make_unique<DifferentialProblemSolver>(params, n);
//         cout << "Entered constructor" << endl;
//         auto a_solver = std::make_unique<DifferentialProblemSolver>(params, n);
//     };

//     SolutionType solve(double tol, int max_iter);
    
//     private:
//     std::unique_ptr<BaseSolver> solver;
//     int m, n, d;
// };

#endif // AUGMENTEDLAGRANGIAN_HPP

// template <class... Args>
// std::unique_ptr<BaseSolver> solverFactory(const std::string & id, Args&&... args) {
//     if (id == "Quadratic"){
//         return std::make_unique<QuadraticSolver>(std::forward<Args>(args)...);
//         }
//     return std::unique_ptr<BaseSolver>(nullptr);
// }
