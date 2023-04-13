//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 10, 1, 1.0};
//+
Point(4) = {0, 10, 1, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 4};
//+
Line(3) = {2, 3};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 3, 2};
//+
Surface(1) = {1};
//+
Transfinite Curve {4, 1, 3, 2} = 10 Using Progression 1;
//+
Transfinite Surface {1};
//+
Physical Curve("left") = {1};
//+
Physical Point("left") += {1, 2};
//+
Physical Point("right") = {4};
