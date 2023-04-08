// Gmsh project created on Wed Apr  5 03:08:16 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {8.660254, 5.0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Physical Point("left", 2) = {1};
//+
Physical Point("right", 3) = {2};
