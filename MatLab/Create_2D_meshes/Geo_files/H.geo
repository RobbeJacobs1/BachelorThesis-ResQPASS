SetFactory("OpenCASCADE");

h = 0.2;    // Size between points

// Define points
Point(1) = {-0.6, -1, 0};
Point(2) = {-0.4, -1, 0};
Point(3) = {-0.4, 1, 0};
Point(4) = {-0.6, 1, 0};

Point(5) = {0.4, -1, 0};
Point(6) = {0.6, -1, 0};
Point(7) = {0.6, 1, 0};
Point(8) = {0.4, 1, 0};

Point(9)  = {-0.4, -0.1, 0};
Point(10) = {0.4, -0.1, 0};
Point(11) = {0.4, 0.1, 0};
Point(12) = {-0.4, 0.1, 0};

// Lines
Line(1) = {1,2};
Line(2) = {2,9};
Line(3) = {9,10};
Line(4) = {10,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,11};
Line(9) = {11,12};
Line(10) = {12,3};
Line(11) = {3,4};
Line(12) = {4,1};

// Surface
Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12};
Plane Surface(1) = {1};

// Global mesh size
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

Mesh.Algorithm = 6;

Mesh 2;
