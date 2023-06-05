// Gmsh project created on Mon Jun  5 02:04:10 2023
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, -5, 10, 10, 0};
//+
Extrude {0, 0, 3} {
  Surface{1}; 
}
//+
Extrude {0, 0, 4} {
  Surface{6}; 
}
//+
Extrude {0, 0, 3} {
  Surface{11}; 
}

//+
Transfinite Curve {23, 20, 28, 27, 19, 25, 17, 11, 3, 9, 2, 1, 7, 4, 12} = 11 Using Progression 1;
//+
Transfinite Curve {26, 24, 22, 21, 5, 6, 8, 10} = 4 Using Progression 1;
//+
Transfinite Curve {13, 14, 16, 18} = 5 Using Progression 1;
//+
Transfinite Volume{1} = {1, 2, 3, 4, 5, 6, 7, 8};
//+
Transfinite Volume{2} = {5, 6, 7, 8, 9, 10, 11, 12};
//+
Transfinite Volume{3} = {9, 10, 11, 12, 13, 14, 15, 16};
//+
Transfinite Surface {2} = {1, 2, 6, 5};
//+
Transfinite Surface {5} = {4, 1, 5, 8};
//+
Transfinite Surface {1} = {1, 2, 3, 4};
//+
Transfinite Surface {3} = {2, 3, 7, 6};
//+
Transfinite Surface {7} = {5, 6, 10, 9};
//+
Transfinite Surface {12} = {9, 10, 14, 13};
//+
Transfinite Surface {8} = {6, 7, 11, 10};
//+
Transfinite Surface {13} = {10, 11, 15, 14};
//+
Transfinite Surface {6} = {5, 6, 7, 8};
//+
Transfinite Surface {11} = {9, 10, 11, 12};
//+
Transfinite Surface {16} = {13, 14, 15, 16};
//+
Transfinite Surface {10} = {8, 5, 9, 12};
//+
Transfinite Surface {15} = {12, 9, 13, 16};
//+
Transfinite Surface {9} = {7, 8, 12, 11};
//+
Transfinite Surface {14} = {11, 12, 16, 15};
//+
Transfinite Curve {23, 15, 7, 1} = 11 Using Progression 1;
//+
Transfinite Curve {11, 3} = 11 Using Progression 1;
//+
Transfinite Surface {4} = {3, 4, 8, 7};
//+
Physical Curve("edge") = {21, 13, 5};
//+
Physical Surface("X0") = {15, 10, 5};
//+
Physical Surface("X1") = {13, 8, 3};
//+
Physical Surface("Y0") = {12, 7, 2};
//+
Physical Surface("Y1") = {14, 9, 4};
//+
Physical Volume("low") = {1};
//+
Physical Volume("mid") = {2};
//+
Physical Volume("top") = {3};
