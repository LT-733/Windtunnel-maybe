# Welcome to my little Wind Tunnel experiment!
Well, kind of a wind tunnel. It really isn't the most accurate thing in the world, because I used the good old [Stam's Stable Fluid Solver](https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf) methodologies to make this tunnel, 
which definitely decreased the accuracy of the project, so treat this more like a visualizer than a scientific feat. The fluids are still descently accurate, though. That is thanks to the actual math I did.

So, this project is essentially a Navier Stokes Solver visualized into a wind tunnel using SDL2. You can import ```.stl``` and ```.obj``` files into the folder and mod the paths of the files read within ```main3.cpp``` to display whatever file you want.
Do notice that the bigger the 3d model, the slower the rendering time. So try using smaller models or your frame rate will easily drop below 5fps.
