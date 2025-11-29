#include <iostream>
#include <fstream>
#include <cmath>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>



using namespace dealii;

// Creating a grid in the domain of a square
void first_grid(){
    // Oggetto che gestisce la griglia (2D), inizialmente vuoto
    Triangulation<2> triangulation;

    // Riempie Triangulation con un ipercubo (in 2D è un quadrato)
    GridGenerator::hyper_cube(triangulation);

    //Raffina la griglia 4 volte
    // 1 -> 4 -> 16 -> 64 -> 256
    triangulation.refine_global(4);

    //Output
    std::ofstream out("grid-1.svg");
    GridOut grid_out;
    grid_out.write_svg(triangulation, out);
    std::cout << "Grid written to grid-1.svg" << std::endl;
}

// Grid Domain
void second_grid(){
    Triangulation<2> triangulation;

    //Crea una "ciambella", con centro, i  due raggi e geometria del guscio in 2D
    const Point<2> center(1,0);
    const double inner_radius = 0.5, outer_radius = 1.0;
    GridGenerator::hyper_shell(triangulation, center, inner_radius, outer_radius, 10);

    //Cicli per Raffinare
    for(unsigned int step = 0; step < 5; ++step){

        //Scorri lungo tutte le celle "attive"
        // Una cella è attiva se non è stata suddivisa in altre celle "figlie"
        for(const auto &cell: triangulation.active_cell_iterators()){

            // Raffiniamo solo le celle che sono all'interno del cerchio preso in considerazione
            for(const auto v : cell -> vertex_indices()){
                const double distance_from_center= center.distance(cell ->vertex(v));

                if(std::fabs(distance_from_center - inner_radius) <= 1e-6 * inner_radius){
                    cell -> set_refine_flag();
                    break;
                }
            }
        }
            // Now we have marked all the cells that we want refined. Now lets refine with the triangulation
        triangulation.execute_coarsening_and_refinement();
    }
    std::ofstream out("grid-2.svg");
    GridOut       grid_out;
    grid_out.write_svg(triangulation, out);
  
    std::cout << "Grid written to grid-2.svg" << std::endl;
}
    



int main(){

    first_grid();
    second_grid();

    return 0;
}