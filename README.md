# Welcome to my little Wind Tunnel experiment!
Well, kind of a wind tunnel. It really isn't the most accurate thing in the world, because I used the good old [Stam's Stable Fluid Solver](https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf) methodologies to make this tunnel, 
which definitely decreased the accuracy of the project, so treat this more like a visualizer than a scientific feat. The fluids are still descently accurate, though. That is thanks to the actual math I did.

So, this project is essentially a Navier Stokes Solver visualized into a wind tunnel using SDL2. You can import ```.stl``` and ```.obj``` files into the folder and mod the paths of the files read within ```main3.cpp``` to display whatever file you want.
Do notice that the bigger the 3d model, the slower the rendering time. So try using smaller models or your frame rate will easily drop below 5fps.

## How to get this project to work
(if you use homebrew you can just ```brew install``` all the required packages)

To get this project to build and run you kind of have to be a developer. You will need LLVM, SDL and defintely openmp. Before all that make sure you have a C++ compiler set up. If you are on Windows (just switch to Linux atp, if you won't I will do it for you so get WSL here) try doing [this](https://learn.microsoft.com/en-us/windows/wsl/install), and if you are on mac, try installing homebrew by doing [this](https://brew.sh/). After getting a C++ compiler, install everything you need and run this: 

```sudo apt update && sudo apt install -y build-essential cmake```

in the root directory and run

```
mkdir build && cd build
cmake ..
cmake --build .
```

And then you can run ```./main3``` and see the aerodynamics of cool 3D models!

<img width="957" height="599" alt="image" src="https://github.com/user-attachments/assets/3272dde8-ba21-4300-8987-a9236be8ccf8" />

### Wait, the model is facing the wrong way!
Yeah... That happens a lot. But the controls are simple. Press I, O, and P to adjust the orientation of your 3D model until the wind is flowing in the right way.

Use + and - keys to increase or decrease the wind speed.

Lastly, drag the image with your mouse to rotate the camera.
