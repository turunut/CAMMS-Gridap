// Gmsh project created on Tue Apr 11 13:03:48 2023
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 50, 10, 10};
//+
Physical Surface("left", 13) = {1};
//+
Physical Point("left", 14) = {2, 1, 3, 4};
//+
Physical Curve("left", 15) = {2, 3, 1, 4};
//+
Physical Curve("right", 16) = {6, 7, 5, 8};
//+
Physical Point("right", 17) = {7, 5, 8, 6};
//+
Physical Surface("right", 18) = {2};
//+
Transfinite Volume{1};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2} Right;
//+
Transfinite Surface {3} Right;
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6} Right;
