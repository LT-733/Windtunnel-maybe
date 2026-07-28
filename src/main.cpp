#include <algorithm>
#include "grid.hpp"
#include "solver2D.hpp"
#include <iostream>
#include <cmath>
#include <SDL.h>

int main(int argc, char* argv[]) {
    if (SDL_Init(SDL_INIT_EVERYTHING) < 0) return -1;

    solver_2d mysolver = solver_2d(200, 100);
    mysolver.set_physical_width(2.0f);
    // Center block obstacle
    for (int y = 40; y < 60; ++y)
    {
        for (int x = 90; x < 110; ++x)
        {
            mysolver.simulation_grid.set_solid(x, y, true);
        }
    }

    SDL_Window* window = SDL_CreateWindow(
        "Wind Tunnel 40 km/h",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        1000,
        500,
        SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE
    );

    SDL_Renderer* renderer = SDL_CreateRenderer(
        window,
        -1,
        SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC
    );

    const float wind_speed_kph = 400.0f;
    const float inlet_speed = wind_speed_kph / 360.0f;
    bool running = true;
    SDL_Event event;

    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT)
                running = false;

            if (event.type == SDL_MOUSEBUTTONDOWN && event.button.clicks == 2) {
                Uint32 flags = SDL_GetWindowFlags(window);
                if (flags & SDL_WINDOW_MAXIMIZED)
                    SDL_RestoreWindow(window);
                else
                    SDL_MaximizeWindow(window);
            }
        }

        constexpr float dt = 1.0f / 60.0f;

        Uint64 start = SDL_GetPerformanceCounter();
        mysolver.step(dt, inlet_speed);
        auto vel = mysolver.simulation_grid.get_velocity(100, 50);
        std::cout << vel.u << " " << vel.v << '\n';
        Uint64 end = SDL_GetPerformanceCounter();

        float ms = (float)(end - start) * 1000.0f /
                   SDL_GetPerformanceFrequency();

        std::cout << "solver ms: " << ms << "\n";

        int win_w, win_h;
        SDL_GetRendererOutputSize(renderer, &win_w, &win_h);

        float sx = win_w / 200.0f;
        float sy = win_h / 100.0f;

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        SDL_SetRenderDrawColor(renderer, 80, 0, 80, 255);

        for (int y = 0; y < 100; ++y) {
            for (int x = 0; x < 200; ++x) {
                if (mysolver.simulation_grid.solid(x, y)) {
                    SDL_Rect r{
                        static_cast<int>(x * sx),
                        static_cast<int>(y * sy),
                        static_cast<int>(std::ceil(sx)),
                        static_cast<int>(std::ceil(sy))
                    };
                    SDL_RenderFillRect(renderer, &r);
                }
            }
        }

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);


        for (int y = 2; y < 98; y += 5) {
            for (int x = 2; x < 198; x += 5) {
                auto vel = mysolver.simulation_grid.get_velocity(x, y);
                float magnitude = std::sqrt(vel.u*vel.u + vel.v*vel.v);
                if (magnitude < 1e-6f)
                    continue;
                float max_speed = inlet_speed * 2.0f; // tune this
                float t = std::clamp(magnitude / max_speed, 0.0f, 1.0f);

                // blue -> cyan -> yellow -> red
                Uint8 r = static_cast<Uint8>(255.0f * t);
                Uint8 g = static_cast<Uint8>(255.0f * (1.0f - std::abs(2.0f*t - 1.0f)));
                Uint8 b = static_cast<Uint8>(255.0f * (1.0f - t));

                SDL_SetRenderDrawColor(renderer, r, g, b, 255);

                float scale = 40.0f;

                int x1 = static_cast<int>(x * sx);
                int y1 = static_cast<int>(y * sy);
                int x2 = static_cast<int>(x1 + vel.u * scale);
                int y2 = static_cast<int>(y1 - vel.v * scale);

                SDL_RenderDrawLine(renderer, x1, y1, x2, y2);

                float angle = std::atan2(
                    static_cast<float>(y1 - y2),
                    static_cast<float>(x2 - x1)
                );

                float size = 6.0f;

                SDL_RenderDrawLine(
                    renderer, x2, y2,
                    x2 - size * std::cos(angle - 0.5f),
                    y2 + size * std::sin(angle - 0.5f)
                );

                SDL_RenderDrawLine(
                    renderer, x2, y2,
                    x2 - size * std::cos(angle + 0.5f),
                    y2 + size * std::sin(angle + 0.5f)
                );
            }
        }
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
