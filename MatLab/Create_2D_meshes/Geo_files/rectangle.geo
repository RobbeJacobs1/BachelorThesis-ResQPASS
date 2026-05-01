SetFactory("OpenCASCADE");

h = 0.2;    // Size between points
Point(1) = {-3, -0.2, 0, 0.05};
Point(2) = {3, -0.2, 0, 0.05};
Point(3) = {3, 0.2, 0, 0.05};
Point(4) = {-3, 0.2, 0, 0.05};

// Create lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Create surface
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Generate mesh
Mesh 2;
