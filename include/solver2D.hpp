#ifndef solver2D
#define solver2D
#include <algorithm>
#include <iostream>
#include <memory>
#include <cmath>
#include <vector>
#include "grid.hpp"

class solver_2d {
public:
    enum Boundary_type {
        scalar = 0,
        horizontal = 1,
        vertical = 2
    };

    solver_2d(int xsize, int ysize)
        : simulation_grid(xsize, ysize, 0),
          viscosity(0.0f),
          inverse_cell_size(1.0f){
    }
    void step(float dt, float inlet_speed){
        simulation_grid.u_prev = simulation_grid.u;
        diffuse(horizontal, simulation_grid.u, simulation_grid.u_prev, viscosity, dt);
        simulation_grid.v_prev = simulation_grid.v;
        diffuse(vertical, simulation_grid.v, simulation_grid.v_prev, viscosity, dt);
        set_inlet_velocity(inlet_speed);
        project();
        simulation_grid.u_prev = simulation_grid.u;
        simulation_grid.v_prev = simulation_grid.v;
        advect(dt);
        project();
    }

    void set_physical_width(float width) {
        if (width > 0.0f) {
            inverse_cell_size = static_cast<float>(simulation_grid.get_xsize()) / width;
        }
    }
    void set_inlet_velocity(float u_val) {
        int ysize = simulation_grid.get_ysize();
        int xsize = simulation_grid.get_xsize();
        for (int y = 1; y < ysize - 1; ++y) {
            if (!simulation_grid.solid(1, y)) {
                int idx1 = 1 + y * xsize;
                int idx0 = 0 + y * xsize;
                simulation_grid.u[1 + y * xsize] = u_val;
                simulation_grid.u[0 + y * xsize] = u_val;
                simulation_grid.u_prev[idx1] = u_val;
                simulation_grid.u_prev[idx0] = u_val;
            }
        }
    }
    Grid<2> simulation_grid;
private:
    float viscosity;
    float inverse_cell_size;
    // @brief: Resolves viscous resistance using implicit relaxation formulas.
    void diffuse(enum Boundary_type type, std::vector<float>& current, const std::vector<float>& prev, float coefficient, float dt) {
        float alpha = dt * coefficient * pow(inverse_cell_size, 2);
        float beta = 1.0f + 4.0f * alpha;
        const int xsize = simulation_grid.get_xsize();
        const int ysize = simulation_grid.get_ysize();
        for (int i = 0; i < 20; ++i) {
            for (int y = 1; y < ysize - 1; ++y) {
                for (int x = 1; x < xsize - 1; ++x) {
                    if (((x + y) % 2 == 0) && !simulation_grid.solid(x, y)) {
                        int idx = x + y * xsize;
                        current[idx] = (prev[idx] + alpha * (current[idx + 1] + current[idx - 1] + current[idx - xsize] + current[idx + xsize])) / beta;
                    }
                }
            }
            set_boundary(type, current);
            for (int y = 1; y < ysize - 1; ++y) {
                for (int x = 1; x < xsize - 1; ++x) {
                    if (((x + y) % 2 != 0) && !simulation_grid.solid(x, y)) {
                        int idx = x + y * xsize;
                        current[idx] = (prev[idx] + alpha * (current[idx + 1] + current[idx - 1] + current[idx - xsize] + current[idx + xsize])) / beta;
                    }
                }
            }
            set_boundary(type, current);
        }
    }
    // @brief: Moves values along the velocity field using backward particle tracing.
    void advect(float dt) {
        int xsize = simulation_grid.get_xsize();
        int ysize = simulation_grid.get_ysize();
        float scaled_dt = dt * inverse_cell_size;

        auto& xcurrent = simulation_grid.u;
        auto& xprev    = simulation_grid.u_prev;
        auto& ycurrent = simulation_grid.v;
        auto& yprev    = simulation_grid.v_prev;

        for (int y = 1; y < ysize - 1; ++y) {
            for (int x = 1; x < xsize - 1; ++x) {
                if (!simulation_grid.solid(x, y)) {
                    int current_idx = x + y * xsize;

                    float u = xprev[current_idx];
                    float v = yprev[current_idx];
                    float xsrc = x - scaled_dt * u;
                    float ysrc = y - scaled_dt * v;

                    float xclamped = std::max(0.5f, std::min(xsize - 1.5f, xsrc));
                    float yclamped = std::max(0.5f, std::min(ysize - 1.5f, ysrc));

                    int x0 = static_cast<int>(std::floor(xclamped));
                    int y0 = static_cast<int>(std::floor(yclamped));

                    float s1 = xclamped - static_cast<float>(x0);
                    float s0 = 1.0f - s1;
                    float t1 = yclamped - static_cast<float>(y0);
                    float t0 = 1.0f - t1;

                    int idx00 = x0 + y0 * xsize;
                    int idx10 = (x0 + 1) + y0 * xsize;
                    int idx01 = x0 + (y0 + 1) * xsize;
                    int idx11 = (x0 + 1) + (y0 + 1) * xsize;

                    xcurrent[current_idx] = t0 * (s0 * xprev[idx00] + s1 * xprev[idx10]) + t1 * (s0 * xprev[idx01] + s1 * xprev[idx11]);

                    ycurrent[current_idx] = t0 * (s0 * yprev[idx00] + s1 * yprev[idx10]) + t1 * (s0 * yprev[idx01] + s1 * yprev[idx11]);
                }
            }
        }
        set_boundary(horizontal, xcurrent);
        set_boundary(vertical, ycurrent);
    }



    // @brief: Enforces mass conservation by removing velocity divergence.
    void project() {
        const int xsize = simulation_grid.get_xsize();
        const int ysize = simulation_grid.get_ysize();

        auto& u = simulation_grid.u;
        auto& v = simulation_grid.v;
        auto& p = simulation_grid.pressure;

        std::vector<float> div(xsize * ysize, 0.0f);

        for (int y = 1; y < ysize - 1; ++y) {
            for (int x = 1; x < xsize - 1; ++x) {
                int idx = x + y * xsize;
                if (!simulation_grid.solid(x, y)){
                    div[idx] = 0.5f *(
                        u[idx + 1] - u[idx - 1] +
                        v[idx + xsize] - v[idx - xsize]
                    );
                } else {
                    div[idx] = 0.0f;
                }
                p[idx] = 0.0f;
            }
        }
        float max_div = 0.0f;

        for (float d : div)
            max_div = std::max(max_div, std::abs(d));

        std::cout << "Max divergence RHS = " << max_div << '\n';

        set_boundary(scalar, div);
        set_boundary(scalar, p);

        constexpr int pressure_iterations = 30;

        for (int k = 0; k < pressure_iterations; ++k){
            // Red cells
            for (int y = 1; y < ysize - 1; ++y){
                for (int x = 1; x < xsize - 1; ++x){
                    if (((x + y) & 1) == 0 && !simulation_grid.solid(x, y)){
                        int idx = x + y * xsize;
                        p[idx] = (
                            -div[idx] +
                            p[idx - 1] +
                            p[idx + 1] +
                            p[idx - xsize] +
                            p[idx + xsize]
                        ) * 0.25f;
                    }
                }
            }
            set_boundary(scalar, p);
            // Black cells
            for (int y = 1; y < ysize - 1; ++y){
                for (int x = 1; x < xsize - 1; ++x){
                    if (((x + y) & 1) == 1 && !simulation_grid.solid(x, y)){
                        int idx = x + y * xsize;

                        p[idx] = (
                            -div[idx] +
                            p[idx - 1] +
                            p[idx + 1] +
                            p[idx - xsize] +
                            p[idx + xsize]
                        ) * 0.25f;
                    }
                }
            }
            set_boundary(scalar, p);
        }
        float residual = 0;

        for (int y=1; y<ysize-1; y++)
            for (int x=1; x<xsize-1; x++){
                int idx = x + y*xsize;

                float lap =
                    p[idx-1]+p[idx+1]+
                    p[idx-xsize]+p[idx+xsize]
                    -4*p[idx];

                residual = std::max(
                    residual,
                    std::abs(lap + div[idx])
                );
            }

        std::cout << "Poisson residual = " << residual << '\n';
        for (int y = 1; y < ysize - 1; ++y) {
            for (int x = 2; x < xsize - 1; ++x) {
                if (!simulation_grid.solid(x, y)){
                    int idx = x + y * xsize;
                    u[idx] -= 0.5f * (p[idx + 1] - p[idx - 1]);
                    v[idx] -= 0.5f * (p[idx + xsize] - p[idx - xsize]);
                }
            }
        }

        set_boundary(horizontal, u);
        set_boundary(vertical, v);
    }
    // @brief: Enforces reflection, inversion, or continuity conditions on the outer perimeter padding layer.
    void set_boundary(enum Boundary_type type, std::vector<float>& field){
        const int xsize = simulation_grid.get_xsize();
        const int ysize = simulation_grid.get_ysize();

        auto idx = [xsize](int x, int y){
            return x + y * xsize;
        };

        for (int y = 1; y < ysize - 1; ++y){
            field[idx(0,y)] = field[idx(1,y)];
            field[idx(xsize-1,y)] = field[idx(xsize-2,y)];
        }

        for (int x = 1; x < xsize - 1; ++x){
            if (type == horizontal or type == vertical){
                field[idx(x,0)] = 0.0f;
                field[idx(x,ysize-1)] = 0.0f;
            }else{
                field[idx(x,0)] = field[idx(x,1)];
                field[idx(x,ysize-1)] = field[idx(x,ysize-2)];
            }
        }

        field[idx(0,0)] = 0.5f * (field[idx(1,0)] + field[idx(0,1)]);
        field[idx(xsize-1,0)] = 0.5f * (field[idx(xsize-2,0)] + field[idx(xsize-1,1)]);
        field[idx(0,ysize-1)] = 0.5f * (field[idx(1,ysize-1)] + field[idx(0,ysize-2)]);
        field[idx(xsize-1,ysize-1)] = 0.5f * (field[idx(xsize-2,ysize-1)] + field[idx(xsize-1,ysize-2)]);
    }
};
#endif