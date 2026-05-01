SetFactory("OpenCASCADE");

h = 0.2;    // Size between points
R = 1;

Disk(1) = {0,0,0,R,R};

Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

Mesh.Algorithm = 6;
Mesh 2;
