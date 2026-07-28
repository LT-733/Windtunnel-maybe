#ifndef solver3D
#define solver3D
#include <algorithm>
#include <iostream>
#include <memory>
#include <cmath>
#include <vector>
#include <chrono>
#include "grid.hpp"
#include "omp.h"

class solver_3d {
public:
    enum Boundary_type {
        scalar = 0,
        horizontal = 1,
        vertical = 2,
        depth = 3
    };

    solver_3d(int xsize, int ysize, int zsize)
        : simulation_grid(xsize, ysize, zsize),
          viscosity(0.0f),
          inlet_velocity(0.0f), inverse_cell_size(1.0f) {
    }
    void step(float dt, float inlet_speed){
        simulation_grid.u_prev = simulation_grid.u;
        // auto timer1 = std::chrono::steady_clock::now();
        diffuse(horizontal, simulation_grid.u, simulation_grid.u_prev, viscosity, dt);
        // auto timer2 = std::chrono::steady_clock::now();
        // std::chrono::duration<float> d = std::chrono::duration_cast<std::chrono::duration<float>>(timer2 - timer1);
        // std::cout<<"diffusing took: "<<d<<" seconds\n";
        simulation_grid.v_prev = simulation_grid.v;
        diffuse(vertical, simulation_grid.v, simulation_grid.v_prev, viscosity, dt);
        simulation_grid.w_prev = simulation_grid.w;
        diffuse(depth, simulation_grid.w, simulation_grid.w_prev, viscosity, dt);
        set_inlet_velocity(inlet_speed);
        // timer1 = std::chrono::steady_clock::now();
        project();
        // timer2 = std::chrono::steady_clock::now();
        // d = std::chrono::duration_cast<std::chrono::duration<float>>(timer2 - timer1);
        // std::cout<<"projection took: "<<d<<" seconds\n";
        simulation_grid.u_prev = simulation_grid.u;
        simulation_grid.v_prev = simulation_grid.v;
        simulation_grid.w_prev = simulation_grid.w;
        // timer1 = std::chrono::steady_clock::now();
        advect(dt);
        // timer2 = std::chrono::steady_clock::now();
        // d = std::chrono::duration_cast<std::chrono::duration<float>>(timer2 - timer1);
        // std::cout<<"advection took: "<<d<<" seconds\n";
        project();
    }

    void set_physical_width(float width) {
        if (width > 0.0f) {
            inverse_cell_size = static_cast<float>(simulation_grid.get_xsize()) / width;
        }
    }
    void set_inlet_velocity(float u_val) {
        inlet_velocity = u_val;
        int ysize = simulation_grid.get_ysize();
        int xsize = simulation_grid.get_xsize();
        int zsize = simulation_grid.get_zsize();
        for (int y = 1; y < ysize - 1; ++y) {
            for (int z = 1; z < zsize-1; ++z) {
                if (!simulation_grid.solid(1, y, z)) {
                    int idx1 = 1 + y * xsize + z * ysize * xsize;
                    int idx0 = 0 + y * xsize + z * ysize * xsize;
                    simulation_grid.u[idx1] = u_val;
                    simulation_grid.u[idx0] = u_val;
                    simulation_grid.u_prev[idx1] = u_val;
                    simulation_grid.u_prev[idx0] = u_val;
                }
            }
        }
    }
    Grid<3> simulation_grid;
private:
    float viscosity;
    float inlet_velocity;
    float inverse_cell_size;
    // @brief: Resolves viscous resistance using implicit relaxation formulas.
    void diffuse(enum Boundary_type type, std::vector<float>& current, const std::vector<float>& prev, float coefficient, float dt) {
        float alpha = dt * coefficient * pow(inverse_cell_size, 2);
        const int xsize = simulation_grid.get_xsize();
        const int ysize = simulation_grid.get_ysize();
        const int zsize = simulation_grid.get_zsize();
        omp_set_num_threads(8);
        #pragma omp parallel
        {
            for (int i = 0; i < 20; ++i) {
                #pragma omp for
                for (int z = 1; z < zsize - 1; ++z) {
                    for (int y = 1; y < ysize - 1; ++y) {
                        for (int x = 1; x < xsize - 1; ++x) {
                            if (((x + y + z) % 2 == 0) && !simulation_grid.solid(x, y, z)) {
                                const int idx = x + y * xsize + z * xsize * ysize;
                                float sum = 0.0f;
                                int count = 0;

                                if (!simulation_grid.solid(x + 1, y, z)) { sum += current[idx + 1]; count++; }
                                if (!simulation_grid.solid(x - 1, y, z)) { sum += current[idx - 1]; count++; }
                                if (!simulation_grid.solid(x, y + 1, z)) { sum += current[idx + xsize]; count++; }
                                if (!simulation_grid.solid(x, y - 1, z)) { sum += current[idx - xsize]; count++; }
                                if (!simulation_grid.solid(x, y, z + 1)) { sum += current[idx + xsize * ysize]; count++; }
                                if (!simulation_grid.solid(x, y, z - 1)) { sum += current[idx - xsize * ysize]; count++; }

                                current[idx] = (prev[idx] + alpha * sum) /
                                               (1.0f + alpha * static_cast<float>(count));                            }
                        }
                    }
                }
                #pragma omp single
                {
                    set_boundary(type, current);
                }
                #pragma omp for
                for (int z = 1; z < zsize - 1; ++z) {
                    for (int y = 1; y < ysize - 1; ++y) {
                        for (int x = 1; x < xsize - 1; ++x) {
                            if (((x + y + z) % 2 != 0) && !simulation_grid.solid(x, y, z)) {
                                const int idx = x + y * xsize + z * xsize * ysize;
                                float sum = 0.0f;
                                int count = 0;

                                if (!simulation_grid.solid(x + 1, y, z)) { sum += current[idx + 1]; count++; }
                                if (!simulation_grid.solid(x - 1, y, z)) { sum += current[idx - 1]; count++; }
                                if (!simulation_grid.solid(x, y + 1, z)) { sum += current[idx + xsize]; count++; }
                                if (!simulation_grid.solid(x, y - 1, z)) { sum += current[idx - xsize]; count++; }
                                if (!simulation_grid.solid(x, y, z + 1)) { sum += current[idx + xsize * ysize]; count++; }
                                if (!simulation_grid.solid(x, y, z - 1)) { sum += current[idx - xsize * ysize]; count++; }

                                current[idx] = (prev[idx] + alpha * sum) /
                                               (1.0f + alpha * static_cast<float>(count));
                            }
                        }
                    }
                }
                #pragma omp single
                {
                    set_boundary(type, current);
                }
            }
        }
    }
    // @brief: Moves values along the velocity field using backward particle tracing.
    void advect(float dt) {
        int xsize = simulation_grid.get_xsize();
        int ysize = simulation_grid.get_ysize();
        int zsize = simulation_grid.get_zsize();
        float scaled_dt = dt * inverse_cell_size;

        auto& xcurrent = simulation_grid.u;
        auto& xprev    = simulation_grid.u_prev;
        auto& ycurrent = simulation_grid.v;
        auto& yprev    = simulation_grid.v_prev;
        auto& zcurrent = simulation_grid.w;
        auto& zprev = simulation_grid.w_prev;
        #pragma omp parallel for collapse(3)
        for (int z = 1; z < zsize - 1; ++z) {
            for (int y = 1; y < ysize - 1; ++y) {
                for (int x = 1; x < xsize - 1; ++x) {
                    if (!simulation_grid.solid(x, y, z)) {
                        int current_idx = x + y * xsize + z * ysize * xsize;

                        float u = xprev[current_idx];
                        float v = yprev[current_idx];
                        float w = zprev[current_idx];
                        float xsrc = x - scaled_dt * u;
                        float ysrc = y - scaled_dt * v;
                        float zsrc = z - scaled_dt * w;

                        float xclamped = std::max(0.5f, std::min(xsize - 1.5f, xsrc));
                        float yclamped = std::max(0.5f, std::min(ysize - 1.5f, ysrc));
                        float zclamped = std::max(0.5f, std::min(zsize - 1.5f, zsrc));

                        int x0 = static_cast<int>(std::floor(xclamped));
                        int y0 = static_cast<int>(std::floor(yclamped));
                        int z0 = static_cast<int>(std::floor(zclamped));

                        float s1 = xclamped - static_cast<float>(x0);
                        float s0 = 1.0f - s1;
                        float t1 = yclamped - static_cast<float>(y0);
                        float t0 = 1.0f - t1;
                        float r1 = zclamped - static_cast<float>(z0);
                        float r0 = 1.0f - r1;
                        if (simulation_grid.solid(x0, y0, z0)) {
                            x0 = x;
                            y0 = y;
                            z0 = z;
                            s1 = t1 = r1 = 0.0f;
                            s0 = t0 = r0 = 1.0f;
                        }

                        const int row_stride = xsize;
                        const int slice_stride = xsize * ysize;

                        int base_idx = x0 + y0 * row_stride + z0 * slice_stride;

                        const int idx000 = base_idx;
                        const int idx100 = simulation_grid.solid(x0 + 1, y0, z0) ? base_idx : (base_idx + 1);
                        const int idx010 = simulation_grid.solid(x0, y0 + 1, z0) ? base_idx : (base_idx + row_stride);
                        const int idx110 = simulation_grid.solid(x0 + 1, y0 + 1, z0) ? base_idx : (base_idx + row_stride + 1);
                        const int idx001 = simulation_grid.solid(x0, y0, z0 + 1) ? base_idx : (base_idx + slice_stride);
                        const int idx101 = simulation_grid.solid(x0 + 1, y0, z0 + 1) ? base_idx : (base_idx + slice_stride + 1);
                        const int idx011 = simulation_grid.solid(x0, y0 + 1, z0 + 1) ? base_idx : (base_idx + slice_stride + row_stride);
                        const int idx111 = simulation_grid.solid(x0 + 1, y0 + 1, z0 + 1) ? base_idx : (base_idx + slice_stride + row_stride + 1);
                        const float x_front = t0 * (s0 * xprev[idx000] + s1 * xprev[idx100]) + t1 * (s0 * xprev[idx010] + s1 * xprev[idx110]);
                        const float y_front = t0 * (s0 * yprev[idx000] + s1 * yprev[idx100]) + t1 * (s0 * yprev[idx010] + s1 * yprev[idx110]);
                        const float z_front = t0 * (s0 * zprev[idx000] + s1 * zprev[idx100]) + t1 * (s0 * zprev[idx010] + s1 * zprev[idx110]);
                        const float x_back  = t0 * (s0 * xprev[idx001] + s1 * xprev[idx101]) + t1 * (s0 * xprev[idx011] + s1 * xprev[idx111]);
                        const float y_back  = t0 * (s0 * yprev[idx001] + s1 * yprev[idx101]) + t1 * (s0 * yprev[idx011] + s1 * yprev[idx111]);
                        const float z_back  = t0 * (s0 * zprev[idx001] + s1 * zprev[idx101]) + t1 * (s0 * zprev[idx011] + s1 * zprev[idx111]);

                        xcurrent[current_idx] = r0 * x_front + r1 * x_back;
                        ycurrent[current_idx] = r0 * y_front + r1 * y_back;
                        zcurrent[current_idx] = r0 * z_front + r1 * z_back;
                    }
                }
            }
        }
        set_boundary(horizontal, xcurrent);
        set_boundary(vertical, ycurrent);
        set_boundary(depth, zcurrent);
    }



    // @brief: Enforces mass conservation by removing velocity divergence.
    void project() {
        const int xsize = simulation_grid.get_xsize();
        const int ysize = simulation_grid.get_ysize();
        const int zsize = simulation_grid.get_zsize();

        auto& u = simulation_grid.u;
        auto& v = simulation_grid.v;
        auto& w = simulation_grid.w;
        auto& p = simulation_grid.pressure;

        std::vector<float> div(xsize * ysize * zsize, 0.0f);

        // auto timer1 = std::chrono::steady_clock::now();
        #pragma omp parallel for collapse(3)
        for (int z = 1; z < zsize - 1; ++z) {
            for (int y = 1; y < ysize - 1; ++y) {
                for (int x = 1; x < xsize - 1; ++x) {
                    int idx = x + y * xsize + z * xsize * ysize;
                    if (!simulation_grid.solid(x, y, z)){
                        const float du = (!simulation_grid.solid(x + 1, y, z) ? u[idx + 1] : 0.0f) - (!simulation_grid.solid(x - 1, y, z) ? u[idx - 1] : 0.0f);
                        const float dv = (!simulation_grid.solid(x, y + 1, z) ? v[idx + xsize] : 0.0f) - (!simulation_grid.solid(x, y - 1, z) ? v[idx - xsize] : 0.0f);
                        const float dw = (!simulation_grid.solid(x, y, z + 1) ? w[idx + xsize * ysize] : 0.0f) - (!simulation_grid.solid(x, y, z - 1) ? w[idx - xsize * ysize] : 0.0f);

                        div[idx] = 0.5f * (du + dv + dw);
                    } else {
                        div[idx] = 0.0f;
                    }
                    p[idx] = 0.0f;
                }
            }
        }
        float max_div = 0.0f;
        // auto timer2 = std::chrono::steady_clock::now();
        // std::chrono::duration<float> d = std::chrono::duration_cast<std::chrono::duration<float>>(timer2 - timer1);
        // std::cout<<"div calculation took: "<<d<<" seconds\n";

        for (float d : div)
            max_div = std::max(max_div, std::abs(d));

        // std::cout << "Max divergence RHS = " << max_div << '\n';

        set_boundary(scalar, div);
        set_boundary(scalar, p);

        constexpr int pressure_iterations = 20;
        // timer1 = std::chrono::steady_clock::now();
        #pragma omp parallel
        {
            for (int k = 0; k < pressure_iterations; ++k){
                // Red cells
                #pragma omp for
                for (int z = 1; z < zsize - 1; ++z) {
                    for (int y = 1; y < ysize - 1; ++y){
                        for (int x = 1; x < xsize - 1; ++x){
                            if (((x + y + z) & 1) == 0 && !simulation_grid.solid(x, y, z)){
                                const int idx = x + y * xsize + z * xsize * ysize;
                                float sum = 0.0f;
                                int count = 0;

                                if (!simulation_grid.solid(x - 1, y, z)) { sum += p[idx - 1]; count++; }
                                if (!simulation_grid.solid(x + 1, y, z)) { sum += p[idx + 1]; count++; }
                                if (!simulation_grid.solid(x, y - 1, z)) { sum += p[idx - xsize]; count++; }
                                if (!simulation_grid.solid(x, y + 1, z)) { sum += p[idx + xsize]; count++; }
                                if (!simulation_grid.solid(x, y, z - 1)) { sum += p[idx - xsize * ysize]; count++; }
                                if (!simulation_grid.solid(x, y, z + 1)) { sum += p[idx + xsize * ysize]; count++; }

                                if (count) p[idx] = (sum - div[idx]) / static_cast<float>(count);
                            }
                        }
                    }
                }
                #pragma omp single
                {
                    set_boundary(scalar, p);
                }
                // Black cells
                #pragma omp for
                for (int z = 1; z < zsize - 1; ++z) {
                    for (int y = 1; y < ysize - 1; ++y){
                        for (int x = 1; x < xsize - 1; ++x){
                            if (((x + y + z) & 1) == 1 && !simulation_grid.solid(x, y, z)){
                                const int idx = x + y * xsize + z * xsize * ysize;
                                float sum = 0.0f;
                                int count = 0;

                                if (!simulation_grid.solid(x - 1, y, z)) { sum += p[idx - 1]; count++; }
                                if (!simulation_grid.solid(x + 1, y, z)) { sum += p[idx + 1]; count++; }
                                if (!simulation_grid.solid(x, y - 1, z)) { sum += p[idx - xsize]; count++; }
                                if (!simulation_grid.solid(x, y + 1, z)) { sum += p[idx + xsize]; count++; }
                                if (!simulation_grid.solid(x, y, z - 1)) { sum += p[idx - xsize * ysize]; count++; }
                                if (!simulation_grid.solid(x, y, z + 1)) { sum += p[idx + xsize * ysize]; count++; }

                                if (count) p[idx] = (sum - div[idx]) / static_cast<float>(count);
                            }
                        }
                    }
                }
                #pragma omp single
                {
                    set_boundary(scalar, p);
                }
            }

        }
        // timer2 = std::chrono::steady_clock::now();
        // d = std::chrono::duration_cast<std::chrono::duration<float>>(timer2 - timer1);
        // std::cout<<"red black iterations took: "<<d<<" seconds\n";
        // float residual = 0;
        //
        // for (int z = 1; z < zsize-1; ++z) {
        //     for (int y=1; y<ysize-1; y++)
        //         for (int x=1; x<xsize-1; x++){
        //             if (!simulation_grid.solid(x, y, z)) {
        //                 int idx = x + y * xsize + z * xsize * ysize;
        //                 float sum = 0.0f;
        //                 int count = 0;
        //
        //                 if (!simulation_grid.solid(x - 1, y, z)) { sum += p[idx - 1]; count++; }
        //                 if (!simulation_grid.solid(x + 1, y, z)) { sum += p[idx + 1]; count++; }
        //                 if (!simulation_grid.solid(x, y - 1, z)) { sum += p[idx - xsize]; count++; }
        //                 if (!simulation_grid.solid(x, y + 1, z)) { sum += p[idx + xsize]; count++; }
        //                 if (!simulation_grid.solid(x, y, z - 1)) { sum += p[idx - xsize * ysize]; count++; }
        //                 if (!simulation_grid.solid(x, y, z + 1)) { sum += p[idx + xsize * ysize]; count++; }
        //
        //                 float lap = sum - static_cast<float>(count) * p[idx];
        //                 residual = std::max(
        //                     residual,
        //                     std::abs(lap + div[idx])
        //                 );
        //             }
        //         }
        // }

        // std::cout << "Poisson residual = " << residual << '\n';
        // timer1 = std::chrono::steady_clock::now();
        #pragma omp parallel for collapse(3)
        for (int z = 1; z < zsize-1; ++z) {
            for (int y = 1; y < ysize - 1; ++y) {
                for (int x = 1; x < xsize - 1; ++x) {
                    if (!simulation_grid.solid(x, y, z)){
                        const int idx = x + y * xsize + z * xsize * ysize;
                        const float pxp = simulation_grid.solid(x + 1, y, z) ? p[idx] : p[idx + 1];
                        const float pxm = simulation_grid.solid(x - 1, y, z) ? p[idx] : p[idx - 1];
                        const float pyp = simulation_grid.solid(x, y + 1, z) ? p[idx] : p[idx + xsize];
                        const float pym = simulation_grid.solid(x, y - 1, z) ? p[idx] : p[idx - xsize];
                        const float pzp = simulation_grid.solid(x, y, z + 1) ? p[idx] : p[idx + xsize * ysize];
                        const float pzm = simulation_grid.solid(x, y, z - 1) ? p[idx] : p[idx - xsize * ysize];

                        u[idx] -= 0.5f * (pxp - pxm);
                        v[idx] -= 0.5f * (pyp - pym);
                        w[idx] -= 0.5f * (pzp - pzm);
                    }
                }
            }
        }
        // timer2 = std::chrono::steady_clock::now();
        // d = std::chrono::duration_cast<std::chrono::duration<float>>(timer2 - timer1);
        // std::cout<<"corrections took: "<<d<<" seconds\n";

        // timer1 = std::chrono::steady_clock::now();
        set_boundary(horizontal, u);
        set_boundary(vertical, v);
        set_boundary(depth, w);
        // timer2 = std::chrono::steady_clock::now();
        // d = std::chrono::duration_cast<std::chrono::duration<float>>(timer2 - timer1);
        // std::cout<<"boundary updates took: "<<d<<" seconds\n";
    }
    // @brief: Enforces reflection, inversion, or continuity conditions on the outer perimeter padding layer.
    void set_boundary(enum Boundary_type type, std::vector<float>& field){
        const int xsize = simulation_grid.get_xsize();
        const int ysize = simulation_grid.get_ysize();
        const int zsize = simulation_grid.get_zsize();

        auto idx = [xsize, ysize](int x, int y, int z) {
            return x + y * xsize + z * xsize * ysize;
        };

        // #pragma omp parallel for collapse(3)
        for (int z = 1; z < zsize - 1; ++z) {
            for (int y = 1; y < ysize - 1; ++y) {
                for (int x = 1; x < xsize - 1; ++x) {
                    if (simulation_grid.solid(x, y, z)) {
                        float sum = 0.0f;
                        int count = 0;
                        if (x - 1 > 0 && !simulation_grid.solid(x - 1, y, z)) {
                            sum += (type == horizontal) ? -field[idx(x - 1, y, z)] : field[idx(x - 1, y, z)];
                            count++;
                        }
                        if (x + 1 < xsize - 1 && !simulation_grid.solid(x + 1, y, z)) {
                            sum += (type == horizontal) ? -field[idx(x + 1, y, z)] : field[idx(x + 1, y, z)];
                            count++;
                        }
                        if (y - 1 > 0 && !simulation_grid.solid(x, y - 1, z)) {
                            sum += (type == vertical) ? -field[idx(x, y - 1, z)] : field[idx(x, y - 1, z)];
                            count++;
                        }
                        if (y + 1 < ysize - 1 && !simulation_grid.solid(x, y + 1, z)) {
                            sum += (type == vertical) ? -field[idx(x, y + 1, z)] : field[idx(x, y + 1, z)];
                            count++;
                        }
                        if (z - 1 > 0 && !simulation_grid.solid(x, y, z - 1)) {
                            sum += (type == depth) ? -field[idx(x, y, z - 1)] : field[idx(x, y, z - 1)];
                            count++;
                        }
                        if (z + 1 < zsize - 1 && !simulation_grid.solid(x, y, z + 1)) {
                            sum += (type == depth) ? -field[idx(x, y, z + 1)] : field[idx(x, y, z + 1)];
                            count++;
                        }
                        if (count > 0) {
                            field[idx(x, y, z)] = sum / static_cast<float>(count);
                        }
                    }
                }
            }
        }

        for (int z = 1; z < zsize - 1; ++z) {
            for (int y = 1; y < ysize - 1; ++y) {
                if (type == horizontal) {
                    field[idx(0, y, z)] = inlet_velocity;
                    field[idx(xsize - 1, y, z)] = 0.0f;
                } else {
                    field[idx(0, y, z)] = field[idx(1, y, z)];
                    field[idx(xsize - 1, y, z)] = field[idx(xsize - 2, y, z)];
                }
            }
        }
        for (int z = 1; z < zsize - 1; ++z) {
            for (int x = 1; x < xsize - 1; ++x) {
                if (type == vertical) {
                    field[idx(x, 0, z)] = 0.0f;
                    field[idx(x, ysize - 1, z)] = 0.0f;
                } else {
                    field[idx(x, 0, z)] = field[idx(x, 1, z)];
                    field[idx(x, ysize - 1, z)] = field[idx(x, ysize - 2, z)];
                }
            }
        }
        for (int y = 1; y < ysize - 1; ++y) {
            for (int x = 1; x < xsize - 1; ++x) {
                if (type == depth) {
                    field[idx(x, y, 0)] = 0.0f;
                    field[idx(x, y, zsize - 1)] = 0.0f;
                } else {
                    field[idx(x, y, 0)] = field[idx(x, y, 1)];
                    field[idx(x, y, zsize - 1)] = field[idx(x, y, zsize - 2)];
                }
            }
        }

        for (int x = 1; x < xsize-1; ++x) {
            if (type == depth or type == vertical){
                field[idx(x, 0, 0)] = 0.0f;
                field[idx(x, ysize - 1, 0)] = 0.0f;
                field[idx(x, 0, zsize - 1)] = 0.0f;
                field[idx(x, ysize - 1, zsize - 1)] = 0.0f;
            }else{
                field[idx(x, 0, 0)] = 0.5f * (field[idx(x, 1, 0)] + field[idx(x, 0, 1)]);
                field[idx(x, ysize - 1, 0)] = 0.5f * (field[idx(x, ysize - 2, 0)] + field[idx(x, ysize - 1, 1)]);
                field[idx(x, 0, zsize - 1)] = 0.5f * (field[idx(x, 1, zsize - 1)] + field[idx(x, 0, zsize - 2)]);
                field[idx(x, ysize - 1, zsize - 1)] = 0.5f * (field[idx(x, ysize - 2, zsize - 1)] + field[idx(x, ysize - 1, zsize - 2)]);
            }
        }
        for (int y = 1; y < ysize - 1; ++y) {
            if (type == depth or type == horizontal){
                field[idx(0, y, 0)] = 0.0f;
                field[idx(xsize - 1, y, 0)] = 0.0f;
                field[idx(0, y, zsize - 1)] = 0.0f;
                field[idx(xsize - 1, y, zsize - 1)] = 0.0f;
            }else {
                field[idx(0, y, 0)] = 0.5f * (field[idx(1, y, 0)] + field[idx(0, y, 1)]);
                field[idx(xsize - 1, y, 0)] = 0.5f * (field[idx(xsize - 2, y, 0)] + field[idx(xsize - 1, y, 1)]);
                field[idx(0, y, zsize - 1)] = 0.5f * (field[idx(1, y, zsize - 1)] + field[idx(0, y, zsize - 2)]);
                field[idx(xsize - 1, y, zsize - 1)] = 0.5f * (field[idx(xsize - 2, y, zsize - 1)] + field[idx(xsize - 1, y, zsize - 2)]);
            }
        }
        for (int z = 1; z < zsize - 1; ++z) {
            if (type == horizontal or type == vertical){
                field[idx(0, 0, z)] = 0.0f;
                field[idx(xsize - 1, 0, z)] = 0.0f;
                field[idx(0, ysize - 1, z)] = 0.0f;
                field[idx(xsize - 1, ysize - 1, z)] = 0.0f;
            }else {
                field[idx(0, 0, z)] = 0.5f * (field[idx(1, 0, z)] + field[idx(0, 1, z)]);
                field[idx(xsize - 1, 0, z)] = 0.5f * (field[idx(xsize - 2, 0, z)] + field[idx(xsize - 1, 1, z)]);
                field[idx(0, ysize - 1, z)] = 0.5f * (field[idx(1, ysize - 1, z)] + field[idx(0, ysize - 2, z)]);
                field[idx(xsize - 1, ysize - 1, z)] = 0.5f * (field[idx(xsize - 2, ysize - 1, z)] + field[idx(xsize - 1, ysize - 2, z)]);
            }
        }
        field[idx(0, 0, 0)] = (field[idx(1, 0, 0)] + field[idx(0, 1, 0)] + field[idx(0, 0, 1)]) / 3.0f;
        field[idx(xsize - 1, 0, 0)] = (field[idx(xsize - 2, 0, 0)] + field[idx(xsize - 1, 1, 0)] + field[idx(xsize - 1, 0, 1)]) / 3.0f;
        field[idx(0, ysize - 1, 0)] = (field[idx(1, ysize - 1, 0)] + field[idx(0, ysize - 2, 0)] + field[idx(0, ysize - 1, 1)]) / 3.0f;
        field[idx(xsize - 1, ysize - 1, 0)] = (field[idx(xsize - 2, ysize - 1, 0)] + field[idx(xsize - 1, ysize - 2, 0)] + field[idx(xsize - 1, ysize - 1, 1)]) / 3.0f;
        field[idx(0, 0, zsize - 1)] = (field[idx(1, 0, zsize - 1)] + field[idx(0, 1, zsize - 1)] + field[idx(0, 0, zsize - 2)]) / 3.0f;
        field[idx(xsize - 1, 0, zsize - 1)] = (field[idx(xsize - 2, 0, zsize - 1)] + field[idx(xsize - 1, 1, zsize - 1)] + field[idx(xsize - 1, 0, zsize - 2)]) / 3.0f;
        field[idx(0, ysize - 1, zsize - 1)] = (field[idx(1, ysize - 1, zsize - 1)] + field[idx(0, ysize - 2, zsize - 1)] + field[idx(0, ysize - 1, zsize - 2)]) / 3.0f;
        field[idx(xsize - 1, ysize - 1, zsize - 1)] = (field[idx(xsize - 2, ysize - 1, zsize - 1)] + field[idx(xsize - 1, ysize - 2, zsize - 1)] + field[idx(xsize - 1, ysize - 1, zsize - 2)]) / 3.0f;
    }
};
#endif