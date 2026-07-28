#ifndef grid
#define grid
#include <vector>
#include <iostream>

class solver_2d;
class solver_3d;

// @brief: this is your grid. You will obtain all the speed and pressure and solid forms in it
template <int Dim>
class Grid {
    friend class solver_2d;
    friend class solver_3d;
    private:
        int xsize, ysize, zsize;
        // There's your x, y, z-dimensional velocity, mapped to u, v, w, respectively
        std::vector<float> u;
        std::vector<float> v;
        std::vector<float> w;
        // There's your logged velocities, at the previous time step
        std::vector<float> u_prev;
        std::vector<float> v_prev;
        std::vector<float> w_prev;
        // Your pressure fields and assessing the obstacle parts
        std::vector<float> pressure;
        std::vector<bool> isSolid;

    public:
        struct Vector3D {
            float u, v, w=0.0f;
        };
        int indexer(int x, int y, int z=0) const {
            // if (x < 0) x = 0;
            // if (x >= xsize) x = xsize-1;
            // if (y < 0) y = 0;
            // if (y >= ysize) y = ysize-1;
            int index = x + y * xsize;
            if constexpr (Dim == 3) {
                // if (z < 0) z = 0;
                // if (z >= zsize) z = zsize - 1;
                index += z * xsize * ysize;
            }
            return index;
        }
        Vector3D get_velocity(int x, int y, int z=0) {
            const int idx = indexer(x, y, z);
            const float u_velocity = u.at(idx);
            const float v_velocity = v.at(idx);
            if constexpr(Dim == 3) {
                const float w_velocity = w.at(idx);
                return {u_velocity, v_velocity, w_velocity};
            }
            return {u_velocity, v_velocity};
        }
        void set_velocity(int x, int y, Vector3D value, int z=0) {
            const int idx = indexer(x, y, z);
            u_prev.at(idx) = u.at(idx);
            v_prev.at(idx) = v.at(idx);
            u.at(idx) = value.u;
            v.at(idx) = value.v;
            if constexpr(Dim==3) {
                w_prev.at(idx) = w.at(idx);
                w.at(idx) = value.w;
            }
        }
        float get_pressure(int x, int y, int z=0) const{
            const int idx = indexer(x, y, z);
            const float p = pressure.at(idx);
            return p;
        }
        bool solid(int x, int y, int z=0) {
            const int idx = indexer(x, y, z);
            return isSolid.at(idx);
        }
        void set_solid(int x, int y, int z=0) {
            const int idx = indexer(x, y, z);
            isSolid.at(idx) = true;
        }
        void reset_solid(int x, int y, int z=0) {
            const int idx = indexer(x, y, z);
            if (isSolid.at(idx) == false) std::cout << "Already False!\n";
            isSolid.at(idx) = false;
        }
        int get_xsize() const { return xsize; }
        int get_ysize() const { return ysize; }
        int get_zsize() const { return zsize; }
        Grid(int x, int y, int z) : xsize(x), ysize(y), zsize(z) {
            int totalSize = 0;
            if constexpr(Dim == 2) {
                zsize = 1;
                totalSize = xsize*ysize;
            }else {
                totalSize = xsize*ysize*zsize;
                w.resize(totalSize);
                w_prev.resize(totalSize);
            }
            u.resize(totalSize);
            v.resize(totalSize);
            u_prev.resize(totalSize);
            v_prev.resize(totalSize);
            pressure.resize(totalSize);
            isSolid.resize(totalSize);
        }
};
#endif

