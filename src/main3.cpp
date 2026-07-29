#include <algorithm>
#include <iostream>
#include <cmath>
#include <SDL.h>
#include <sstream>
#include "grid.hpp"
#include "solver3D.hpp"
#include <fstream>


const int WIDTH = 1200;
const int HEIGHT = 800;


struct Vec3
{
    float x;
    float y;
    float z;
};

struct Triangle {
    Vec3 a,b,c;
};

struct Point2D
{
    int x;
    int y;
};


// Projection controls
float CAMERA_DISTANCE = 10.0f;
float SCREEN_SCALE = 500.0f;


// Add these globals near CAMERA_DISTANCE
float camera_pitch = 0.0f; // Up/Down rotation
float camera_yaw = 0.0f;   // Left/Right rotation
float drag_sensitivity = 0.0015f;

std::vector<Triangle> load_stl(const std::string& filename);
std::vector<Triangle> load_obj(const std::string& filename);
void normalize_mesh(std::vector<Triangle>& triangles);

Point2D project(Vec3 p, float scale)
{
    float camera_distance = 10.0f;

    float cosY = std::cos(camera_yaw);
    float sinY = std::sin(camera_yaw);

    float x1 = p.x * cosY - p.z * sinY;
    float z1 = p.x * sinY + p.z * cosY;


    float cosX = std::cos(camera_pitch);
    float sinX = std::sin(camera_pitch);

    float y2 = p.y * cosX - z1 * sinX;
    float z2 = p.y * sinX + z1 * cosX;


    // prevent perspective explosion
    float depth = z2 + camera_distance;

    if(depth < 1.0f)
        depth = 1.0f;


    return {
        WIDTH / 2 + static_cast<int>((x1 / depth) * scale),
        HEIGHT / 2 - static_cast<int>((y2 / depth) * scale)
    };
}


// Converts grid coordinates into a larger world
Vec3 grid_to_world(int x, int y, int z,
                   int xs, int ys, int zs)
{
    float scale = 8.0f;

    return {
        (x - xs / 2.0f) / scale,
        (y - ys / 2.0f) / scale,
        (z - zs / 2.0f) / scale
    };
}

std::vector<Triangle> load_stl(const std::string& filename){
    std::ifstream file(filename, std::ios::binary);

    if (!file)
        throw std::runtime_error("Failed to open " + filename);

    char header[80];
    file.read(header, 80);

    uint32_t triangle_count;
    file.read(reinterpret_cast<char*>(&triangle_count), 4);

    std::vector<Triangle> triangles;
    triangles.reserve(triangle_count);

    for (uint32_t i = 0; i < triangle_count; ++i){
        float normal[3];
        float v[9];
        uint16_t attribute;

        file.read(reinterpret_cast<char*>(normal), sizeof(normal));
        file.read(reinterpret_cast<char*>(v), sizeof(v));
        file.read(reinterpret_cast<char*>(&attribute), sizeof(attribute));

        if (!file)
            throw std::runtime_error("Corrupt STL file.");

        triangles.push_back({
            {v[0], v[1], v[2]},
            {v[3], v[4], v[5]},
            {v[6], v[7], v[8]}
        });
    }

    return triangles;
}

std::vector<Triangle> load_obj(const std::string& filename){
    std::ifstream file(filename);
    if (!file)
        throw std::runtime_error("Failed to open " + filename);
    std::vector<Vec3> vertices;
    std::vector<Triangle> triangles;
    std::string line;
    while (std::getline(file, line)){
        std::stringstream ss(line);

        std::string type;
        ss >> type;

        if (type == "v"){
            Vec3 v;
            ss >> v.x >> v.y >> v.z;
            vertices.push_back(v);
        }else if (type == "f"){
            std::string s1,s2,s3;

            ss >> s1 >> s2 >> s3;

            int a = std::stoi(s1.substr(0,s1.find('/')));
            int b = std::stoi(s2.substr(0,s2.find('/')));
            int c = std::stoi(s3.substr(0,s3.find('/')));
            triangles.push_back({
                vertices[a-1],
                vertices[b-1],
                vertices[c-1]
            });
        }
    }
    return triangles;
}

std::vector<Triangle> load_mesh(const std::string& filename){
    auto extension = filename.substr(filename.find_last_of('.') + 1);
    if (extension == "stl")
        return load_stl(filename);
    if (extension == "obj")
        return load_obj(filename);
    throw std::runtime_error("Unsupported mesh format");
}
float mesh_rot_x = 0.0f;
float mesh_rot_y = 0.0f;
float mesh_rot_z = 0.0f;

void orient_mesh_to_wind(std::vector<Triangle>& triangles) {
    float cx = std::cos(mesh_rot_x), sx = std::sin(mesh_rot_x);
    float cy = std::cos(mesh_rot_y), sy = std::sin(mesh_rot_y);
    float cz = std::cos(mesh_rot_z), sz = std::sin(mesh_rot_z);

    for (auto& t : triangles) {
        auto rot = [&](Vec3 v) {
            float y1 = v.y * cx - v.z * sx;
            float z1 = v.y * sx + v.z * cx;

            float x2 = v.x * cy + z1 * sy;
            float z2 = -v.x * sy + z1 * cy;

            float x3 = x2 * cz - y1 * sz;
            float y3 = x2 * sz + y1 * cz;

            return Vec3{x3, y3, z2};
        };
        t.a = rot(t.a);
        t.b = rot(t.b);
        t.c = rot(t.c);
    }
}
void normalize_mesh(std::vector<Triangle>& triangles,float target_size){
    float minx=FLT_MAX,miny=FLT_MAX,minz=FLT_MAX;
    float maxx=-FLT_MAX,maxy=-FLT_MAX,maxz=-FLT_MAX;

    for(auto& t:triangles){
        for(auto v:{t.a,t.b,t.c}){
            minx=std::min(minx,v.x);
            maxx=std::max(maxx,v.x);
            miny=std::min(miny,v.y);
            maxy=std::max(maxy,v.y);
            minz=std::min(minz,v.z);
            maxz=std::max(maxz,v.z);
        }
    }

    float size=std::max({maxx-minx,maxy-miny,maxz-minz});
    float scale=target_size/size;

    for(auto& t:triangles){
        t.a.x*=scale;t.a.y*=scale;t.a.z*=scale;
        t.b.x*=scale;t.b.y*=scale;t.b.z*=scale;
        t.c.x*=scale;t.c.y*=scale;t.c.z*=scale;
    }
}

int main(int argc, char* argv[])
{
    const std::vector<Triangle> base_triangles = load_mesh("../porsche_gt3rs_spoiler.obj");
    std::vector<Triangle> triangles = base_triangles;
    std::cout << "triangles: "
              << triangles.size()
              << "\n";


    if(triangles.empty()){
        std::cout << "No triangles loaded\n";
        return 1;
    }


    // shrink STL coordinates
    float mesh_scale = 0.50f;

    for(auto& t : triangles){
        t.a.x *= mesh_scale;
        t.a.y *= mesh_scale;
        t.a.z *= mesh_scale;

        t.b.x *= mesh_scale;
        t.b.y *= mesh_scale;
        t.b.z *= mesh_scale;

        t.c.x *= mesh_scale;
        t.c.y *= mesh_scale;
        t.c.z *= mesh_scale;
    }
    normalize_mesh(triangles,6.0f);
    orient_mesh_to_wind(triangles);

    float minx=FLT_MAX,miny=FLT_MAX,minz=FLT_MAX;
    float maxx=-FLT_MAX,maxy=-FLT_MAX,maxz=-FLT_MAX;

    for(auto& t:triangles){
        for(auto v:{t.a,t.b,t.c}){
            minx=std::min(minx,v.x);
            maxx=std::max(maxx,v.x);
            miny=std::min(miny,v.y);
            maxy=std::max(maxy,v.y);
            minz=std::min(minz,v.z);
            maxz=std::max(maxz,v.z);
        }
    }

    float cx=(minx+maxx)*0.5f;
    float cy=(miny+maxy)*0.5f;
    float cz=(minz+maxz)*0.5f;

    for(auto& t:triangles){
        t.a.x-=cx;t.a.y-=cy;t.a.z-=cz;
        t.b.x-=cx;t.b.y-=cy;t.b.z-=cz;
        t.c.x-=cx;t.c.y-=cy;t.c.z-=cz;
    }
    SDL_Init(SDL_INIT_EVERYTHING);


    solver_3d mysolver(
        80,
        50,
        50
    );


    int solid_count = 0;


    auto voxelize = [&](Vec3 v){
        int x = 40 + (int)(v.x * 8.0f);
        int y = 25 + (int)(v.y * 8.0f);
        int z = 25 + (int)(v.z * 8.0f);

        if(x >= 0 && x < 80 &&
           y >= 0 && y < 50 &&
           z >= 0 && z < 50){
            mysolver.simulation_grid.set_solid(x,y,z);
            solid_count++;
        }
    };


    for(auto& t : triangles){
        const int steps = 10;
        for(int i=0;i<=steps;i++){
            for(int j=0;j<=steps-i;j++){
                float a = (float)i / steps;
                float b = (float)j / steps;
                float c = 1.0f - a - b;

                Vec3 p {
                    t.a.x*a +
                    t.b.x*b +
                    t.c.x*c,

                    t.a.y*a +
                    t.b.y*b +
                    t.c.y*c,

                    t.a.z*a +
                    t.b.z*b +
                    t.c.z*c
                };

                voxelize(p);
            }
        }
    }


    std::cout
        << "solid cells: "
        << solid_count
        << "\n";


    mysolver.set_physical_width(2.0f);

    constexpr float dt = 1.0f / 60.0f;
    float scale = 800.0f;
    float inlet_speed = 3.0f;
    float speed_increment = 0.5f;


    SDL_Window* window = SDL_CreateWindow(
        "3D Wind Tunnel",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        WIDTH,
        HEIGHT,
        SDL_WINDOW_SHOWN |
        SDL_WINDOW_RESIZABLE
    );


    SDL_Renderer* renderer =
        SDL_CreateRenderer(
            window,
            -1,
            SDL_RENDERER_ACCELERATED
        );



    bool running = true;
    SDL_Event event;




    while(running){
        while(SDL_PollEvent(&event)){
            if(event.type == SDL_QUIT)
                running = false;

            if(event.type == SDL_MOUSEBUTTONDOWN &&
               event.button.clicks == 2){
                Uint32 flags = SDL_GetWindowFlags(window);

                if(flags & SDL_WINDOW_MAXIMIZED)
                    SDL_RestoreWindow(window);
                else
                    SDL_MaximizeWindow(window);
            }
            if (event.type == SDL_MOUSEWHEEL) {
                scale += static_cast<float>(event.wheel.y);
            }

            if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_LEFT){
                bool dragging = true;

                // Clear out any old relative mouse delta accumulated before clicking
                int flush_x, flush_y;
                SDL_GetRelativeMouseState(&flush_x, &flush_y);

                while (dragging){
                    SDL_PumpEvents(); // Pull latest OS hardware signals

                    // Break out instantly if they let go of left click
                    Uint32 mouse_state = SDL_GetMouseState(nullptr, nullptr);
                    if (!(mouse_state & SDL_BUTTON(SDL_BUTTON_LEFT))) {
                        dragging = false;
                        break;
                    }

                    // Fetch how far the hardware moved since the last iteration
                    int xrel, yrel;
                    SDL_GetRelativeMouseState(&xrel, &yrel);

                    // Map mouse changes to camera rotation angles
                    camera_yaw   += xrel * drag_sensitivity;
                    camera_pitch -= yrel * drag_sensitivity; // Inverted for natural look/feel

                    // Cap pitch so the camera doesn't flip completely upside down
                    camera_pitch = std::clamp(camera_pitch, -1.0f, 1.0f);

                    // Minor delay to keep this thread from spinning at 100% capacity
                    SDL_Delay(1);
                }
            }
            if(event.type == SDL_KEYDOWN)
            {
                if(event.key.keysym.sym == SDLK_EQUALS ||
                   event.key.keysym.sym == SDLK_PLUS)
                {
                    inlet_speed += speed_increment;
                    // std::cout << "Wind speed: "
                    //           << inlet_speed
                    //           << "\n";
                }

                if(event.key.keysym.sym == SDLK_MINUS)
                {
                    inlet_speed -= speed_increment;

                    if(inlet_speed < 0.0f)
                        inlet_speed = 0.0f;

                    // std::cout << "Wind speed: "
                    //           << inlet_speed
                    //           << "\n";
                }
                if (event.key.keysym.sym == SDLK_i) {
                    mesh_rot_x += M_PI / 2.0f;
                }
                if (event.key.keysym.sym == SDLK_o) {
                    mesh_rot_y += M_PI / 2.0f;
                }
                if (event.key.keysym.sym == SDLK_p) {
                    mesh_rot_z += M_PI / 2.0f;
                }

                if (event.key.keysym.sym == SDLK_i || event.key.keysym.sym == SDLK_o || event.key.keysym.sym == SDLK_p) {
                    triangles = base_triangles; // Reset directly from RAM copy

                    for (auto& t : triangles) {
                        t.a.x *= mesh_scale; t.a.y *= mesh_scale; t.a.z *= mesh_scale;
                        t.b.x *= mesh_scale; t.b.y *= mesh_scale; t.b.z *= mesh_scale;
                        t.c.x *= mesh_scale; t.c.y *= mesh_scale; t.c.z *= mesh_scale;
                    }
                    normalize_mesh(triangles, 6.0f);
                    orient_mesh_to_wind(triangles);

                    float minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX;
                    float maxx = -FLT_MAX, maxy = -FLT_MAX, maxz = -FLT_MAX;
                    for (auto& t : triangles) {
                        for (auto v : {t.a, t.b, t.c}) {
                            minx = std::min(minx, v.x); maxx = std::max(maxx, v.x);
                            miny = std::min(miny, v.y); maxy = std::max(maxy, v.y);
                            minz = std::min(minz, v.z); maxz = std::max(maxz, v.z);
                        }
                    }
                    cx = (minx + maxx) * 0.5f; cy = (miny + maxy) * 0.5f; cz = (minz + maxz) * 0.5f;
                    for (auto& t : triangles) {
                        t.a.x -= cx; t.a.y -= cy; t.a.z -= cz;
                        t.b.x -= cx; t.b.y -= cy; t.b.z -= cz;
                        t.c.x -= cx; t.c.y -= cy; t.c.z -= cz;
                    }

                    mysolver = solver_3d(80, 50, 50);
                    solid_count = 0;
                    for (auto& t : triangles) {
                        const int steps = 10;
                        for (int i = 0; i <= steps; i++) {
                            for (int j = 0; j <= steps - i; j++) {
                                float a = (float)i / steps; float b = (float)j / steps; float c = 1.0f - a - b;
                                Vec3 p{
                                    t.a.x * a + t.b.x * b + t.c.x * c,
                                    t.a.y * a + t.b.y * b + t.c.y * c,
                                    t.a.z * a + t.b.z * b + t.c.z * c
                                };
                                voxelize(p);
                            }
                        }
                    }
                }
            }
        }
        Uint64 t0 = SDL_GetPerformanceCounter();
        mysolver.step(dt, inlet_speed);
        std::stringstream title;
        title<<"3D Wind Tunnel, Inlet Speed: "<< inlet_speed;
        SDL_SetWindowTitle(window,title.str().c_str());
        Uint64 t1 = SDL_GetPerformanceCounter();
        SDL_SetRenderDrawColor(
            renderer,
            0,
            0,
            0,
            255
        );

        SDL_RenderClear(renderer);

        SDL_SetRenderDrawColor(renderer,0,255,0,255);

        for(float y=-2;y<=2;y+=0.2f)
        {
            Point2D a = project({-4,y,0},scale);
            Point2D b = project({ 4,y,0},scale);

            SDL_RenderDrawLine(
                renderer,
                a.x,a.y,
                b.x,b.y
            );
        }
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);

        // X axis
        Point2D x1 = project({-3,0,0}, scale);
        Point2D x2 = project({ 3,0,0}, scale);

        SDL_RenderDrawLine(renderer, x1.x, x1.y, x2.x, x2.y);

        int xs = mysolver.simulation_grid.get_xsize();
        int ys = mysolver.simulation_grid.get_ysize();
        int zs = mysolver.simulation_grid.get_zsize();

        for (int z = 0; z < zs; ++z){
            for (int y = 0; y < ys; ++y){
                for (int x = 0; x < xs; ++x){
                    if (!mysolver.simulation_grid.solid(x, y, z))
                        continue;

                    Vec3 c = grid_to_world(x, y, z, xs, ys, zs);

                    float s = 0.06f;   // half cube size

                    Vec3 v[8] ={
                        {c.x-s,c.y-s,c.z-s},
                        {c.x+s,c.y-s,c.z-s},
                        {c.x+s,c.y+s,c.z-s},
                        {c.x-s,c.y+s,c.z-s},

                        {c.x-s,c.y-s,c.z+s},
                        {c.x+s,c.y-s,c.z+s},
                        {c.x+s,c.y+s,c.z+s},
                        {c.x-s,c.y+s,c.z+s}
                    };

                    Point2D p[8];
                    for(int i=0;i<8;i++)
                        p[i]=project(v[i],scale);

                    SDL_SetRenderDrawColor(renderer,0,80,200,255);

                    const int edges[12][2]={
                        {0,1},{1,2},{2,3},{3,0},
                        {4,5},{5,6},{6,7},{7,4},
                        {0,4},{1,5},{2,6},{3,7}
                    };

                    for(auto &e : edges){
                        SDL_RenderDrawLine(
                            renderer,
                            p[e[0]].x,p[e[0]].y,
                            p[e[1]].x,p[e[1]].y
                        );
                    }
                }
            }
        }

        SDL_SetRenderDrawColor(renderer, 255,255,255,255);

        for(auto& t : triangles){
            Point2D a = project(t.a, scale);
            Point2D b = project(t.b, scale);
            Point2D c = project(t.c, scale);

            SDL_RenderDrawLine(renderer,
                a.x,a.y,
                b.x,b.y);

            SDL_RenderDrawLine(renderer,
                b.x,b.y,
                c.x,c.y);

            SDL_RenderDrawLine(renderer,
                c.x,c.y,
                a.x,a.y);
        }

        for(int z=20; z<40; z+=5){
            for(int y=30; y<50; y+=5){
                for(int x=20; x<80; x+=5){

                    if(mysolver.simulation_grid.solid(x,y,z))
                        continue;


                    auto vel =
                        mysolver.simulation_grid.get_velocity(x,y,z);


                    float speed =
                        std::sqrt(
                            vel.u*vel.u +
                            vel.v*vel.v +
                            vel.w*vel.w
                        );


                    if(speed < 0.01f)
                        continue;

                    float max_speed = 1.5f; // tune this
                    float t = std::clamp(speed / max_speed, 0.0f, 1.0f);

                    // blue -> cyan -> yellow -> red
                    Uint8 r = static_cast<Uint8>(255.0f * t);
                    Uint8 g = static_cast<Uint8>(255.0f * (1.0f - std::abs(2.0f*t - 1.0f)));
                    Uint8 blue = static_cast<Uint8>(255.0f * (1.0f - t));

                    SDL_SetRenderDrawColor(
                        renderer,
                        r,
                        g,
                        blue,
                        255
                    );

                    Vec3 start =
                        grid_to_world(
                            x,y,z,
                            xs,ys,zs
                        );


                    Vec3 finish =
                    {
                        start.x + vel.u * 0.3f,
                        start.y + vel.v * 0.3f,
                        start.z + vel.w * 0.3f
                    };


                    Point2D a = project(start, scale);
                    Point2D b = project(finish, scale);


                    SDL_RenderDrawLine(
                        renderer,
                        a.x,
                        a.y,
                        b.x,
                        b.y
                    );


                    float angle =
                        std::atan2(
                            b.y-a.y,
                            b.x-a.x
                        );


                    float arrow = 8.0f;


                    SDL_RenderDrawLine(
                        renderer,
                        b.x,
                        b.y,
                        b.x - arrow*cos(angle-0.5f),
                        b.y - arrow*sin(angle-0.5f)
                    );


                    SDL_RenderDrawLine(
                        renderer,
                        b.x,
                        b.y,
                        b.x - arrow*cos(angle+0.5f),
                        b.y - arrow*sin(angle+0.5f)
                    );
                }
            }
        }


        Uint64 t2 = SDL_GetPerformanceCounter();


        SDL_RenderPresent(renderer);


        Uint64 t3 = SDL_GetPerformanceCounter();


        // float solver_ms =
        //     (float)(t1-t0) *
        //     1000.0f /
        //     SDL_GetPerformanceFrequency();
        //
        // float draw_ms =
        //     (float)(t2-t1) *
        //     1000.0f /
        //     SDL_GetPerformanceFrequency();
        //
        // float present_ms =
        //     (float)(t3-t2) *
        //     1000.0f /
        //     SDL_GetPerformanceFrequency();
        //
        //
        // std::cout
        //     << "solver: "
        //     << solver_ms
        //     << " ms | draw: "
        //     << draw_ms
        //     << " ms | present: "
        //     << present_ms
        //     << " ms\n";
    }



    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();


    return 0;
}