#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
 
#include <deal.II/fe/fe_q.h>
 
#include <deal.II/dofs/dof_tools.h>
 
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
 
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
 
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
 
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
 
// PDE di Poisson su un quadrato [-1, 1]^2
/*
- A mesh on which to define shape functions. You have already seen how to generate and manipulate the objects that describe meshes in step-1 and step-2.
- A finite element that describes the shape functions we want to use on the reference cell (which in deal.II is always the unit interval [0,1], the unit square [0,1]^2 or the unit cube [0,1]^3, depending on which 
space dimension you work in). In step-2, we had already used an object of type FE_Q<2>, which denotes the usual Lagrange elements that define shape functions by interpolation on support points. 
The simplest one is FE_Q<2>(1), which uses polynomial degree 1. In 2d, these are often referred to as bilinear, since they are linear in each of the two coordinates of the reference cell. (In 1d, they would be linear and in 3d tri-linear; 
however, in the deal.II documentation, we will frequently not make this distinction and simply always call these functions "linear".)
- A DoFHandler object that enumerates all the degrees of freedom on the mesh, taking the reference cell description the finite element object provides as the basis. You've also already seen how to do this in step-2.
- A mapping that tells how the shape functions on the real cell are obtained from the shape functions defined by the finite element class on the reference cell. By default, unless you explicitly say otherwise, deal.II 
will use a (bi-, tri-)linear mapping for this, so in most cases you don't have to worry about this step.
*/

using namespace dealii;

class Step3{
    public:
        Step3();
        void run();

    private:
        void make_grid();
        void setup_system();
        void assemble_system();
        void solve();
        void output_results() const;

        Triangulation<2> triangulation;
        FE_Q<2> fe;
        DoFHandler <2> dof_handler;

        SparsityPattern sparsity_pattern;
        SparseMatrix<double> system_matrix;
        Vector<double> solution;
        Vector<double> system_rhs;
};

Step3::Step3()
    : fe(1)
    , dof_handler(triangulation)
    {}

void Step3::make_grid(){
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(5);

    std::cout << "Number of active cells: " << triangulation.n_active_cells() << std::endl;
}

void Step3::setup_system(){
    dof_handler.distribute_dofs(fe);
    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}

void Step3::assemble_system(){
    const QGauss<2> quadrature_formula(fe.degree + 1);

    FEValues<2> fe_values(fe, quadrature_formula, update_values| update_gradients | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Assemblaggio locale
    for(const auto &cell : dof_handler.active_cell_iterators()){
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs = 0;
        for(const unsigned int q_index : fe_values.quadrature_point_indices()){
            for(const unsigned int i : fe_values.dof_indices())
                for(const unsigned int j :fe_values.dof_indices())
                    cell_matrix(i, j) += (fe_values.shape_grad(i, q_index) *
                                         fe_values.shape_grad(j, q_index) *
                                         fe_values.JxW(q_index));
            for(const unsigned int i : fe_values.dof_indices())
                cell_rhs(i) += (fe_values.shape_value(i,q_index)*
                               1.*
                               fe_values.JxW(q_index));                                 
        }
        cell -> get_dof_indices(local_dof_indices);
        for(const unsigned int i : fe_values.dof_indices())
            for(const unsigned int j : fe_values.dof_indices())
                system_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
        
        for(const unsigned int i : fe_values.dof_indices())
            system_rhs(local_dof_indices[i]) += cell_rhs(i);        
    }

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler, types::boundary_id(0), Functions::ZeroFunction<2>(), boundary_values);

    MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs);
}

void Step3::solve(){
    SolverControl solver_control(1000, 1e-6 * system_rhs.l2_norm());
    SolverCG<Vector<double>> solver(solver_control);
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
    
    std::cout << solver_control.last_step() << "CG iterations needed to obtain convergence." << std::endl;
}

void Step3::output_results() const{
    DataOut<2> data_out;
    
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");

    data_out.build_patches();
    
    const std::string filename = "solution.vtk";
    std::ofstream     output(filename);
    data_out.write_vtk(output);
    std::cout << "Output written to " << filename << std::endl;
}

void Step3::run(){
    make_grid();
    setup_system();
    assemble_system();
    solve();
    output_results();
}


int main(){

    Step3 laplace_problem;
    laplace_problem.run();

    return 0;
}