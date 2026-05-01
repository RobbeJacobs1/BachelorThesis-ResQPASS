SetFactory("OpenCASCADE");

h = 0.2;    // Size between points
r = 1.0;          // Radius
xc = 0; yc = 0; zc = 0;   // Center

// Mesh size settings
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

// Create solid sphere
Sphere(1) = {xc, yc, zc, r};

