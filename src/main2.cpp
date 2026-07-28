#include <iostream>
#include <vector>
#include <iomanip>
#include "solver3D.hpp"
#include <chrono>

int main() {
    // 1. Initialize Domain Sizes
    const int XSIZE = 16;
    const int YSIZE = 8;
    const int ZSIZE = 8;

    std::cout << "Initializing solver [" << XSIZE << "x" << YSIZE << "x" << ZSIZE << "]...\n";
    solver_3d solver(XSIZE, YSIZE, ZSIZE);

    // Configure solver physical baseline properties
    solver.set_physical_width(4.0f);

    // 2. Mock Vehicle Solid Geometry using your explicit method
    std::cout << "Injecting flat obstacle profile at Y=1, Z=[3,4], X=[4..9]...\n";
    for (int x = 4; x <= 9; ++x) {
        for (int z = 3; z <= 4; ++z) {
            solver.simulation_grid.set_solid(x, 1, z);
            solver.simulation_grid.set_solid(x, 2, z);
        }
    }

    // 3. Execution Driver Parameters
    const float dt = 0.01f;
    const float inlet_speed = 10.0f;
    const int total_steps = 5;

    std::cout << "\nStarting Simulation Verification Loop (" << total_steps << " steps)...\n";
    std::cout << "========================================================\n";

    for (int step_idx = 1; step_idx <= total_steps; ++step_idx) {
        std::cout << "\n--- TIMESTEP " << step_idx << " ---\n";

        // Execute the single complete fractional step pipeline
        auto start = std::chrono::high_resolution_clock::now();

        solver.step(dt, inlet_speed);

        auto end = std::chrono::high_resolution_clock::now();

        double ms =
            std::chrono::duration<double, std::milli>(end - start).count();

        std::cout << "Solver step time: "
                  << ms
                  << " ms\n";
        // Peek using your actual built-in public getter API
        auto vel = solver.simulation_grid.get_velocity(3, 1, 3);
        std::cout << "Fluid velocity immediately in front of obstacle (3,1,3): "
                  << "U=" << vel.u << ", "
                  << "V=" << vel.v << ", "
                  << "W=" << vel.w << "\n";
    }

    std::cout << "\n========================================================\n";
    std::cout << "Simulation loop run completed cleanly.\n";

    // 4. Print out a slice diagram using your indexer/getter structures
    int target_z = 3;
    int xsize = solver.simulation_grid.get_xsize();
    int ysize = solver.simulation_grid.get_ysize();

    std::cout << "\n--- Velocity U Slice at Z = " << target_z << " ---\n";
    for (int y = ysize - 1; y >= 0; --y) {
        std::cout << "Y=" << std::setw(2) << y << " | ";
        for (int x = 0; x < xsize; ++x) {
            if (solver.simulation_grid.solid(x, y, target_z)) {
                std::cout << " [SS]  ";
            } else {
                auto v_cell = solver.simulation_grid.get_velocity(x, y, target_z);
                std::cout << std::setw(6) << std::fixed << std::setprecision(2) << v_cell.u << " ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "-------------------------------------\n";

    return 0;
}