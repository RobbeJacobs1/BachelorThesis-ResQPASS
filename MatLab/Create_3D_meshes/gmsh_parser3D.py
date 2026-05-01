import numpy as np
from scipy.io import savemat
import subprocess
import itertools
# Note that this code assumes that the Gmsh executable is available in the system's PATH and that the 
# necessary directories for geo files, mesh files, and parsed meshes exist. The code generates 3D meshes 
# from the specified geo files with varying refinement levels, parses the generated mesh files to extract 
# nodes and elements, and saves the parsed data in .mat format for further use.

# Define the list of geo files and their corresponding refinement levels (h values)
filenames = []
for h in np.linspace(0.05, 0.20, 16):
    filenames.append((f"Geo_files/sphere.geo", round(h, 2)))

# Define the directories for the geo files, mesh files, and parsed mesh files
GEO_DIR = "Geo_files/"
MESH_DIR = "Mesh_files/"
PARSED_DIR = "Parsed_meshes/"


def generate_meshes(file: str, refinement: float|None) -> str:
    # This function generates a 3D mesh from a geo file using Gmsh, with an optional refinement level (h value). 
    # It modifies the geo file to set the desired refinement level, runs the Gmsh command to generate the mesh, 
    # and returns the path to the generated mesh file.
    mesh_filename = f"{MESH_DIR}{file[len(GEO_DIR):-4]}_{refinement:.2f}.msh" if refinement is not None else f"{MESH_DIR}{file[len(GEO_DIR):-4]}.msh"
    if refinement is not None:
        # Get the information from the geo file
        with open(file, "r") as f:
            data = f.readlines()

        # Change point_size to the refinement value
        data[2] = f"h = {refinement};    // Size between points\n"

        # Write the modified data back to the geo file
        with open(file, "w") as f:
            f.writelines(data)      
    
    # Run the Gmsh command to generate the mesh
    out = subprocess.run([
        "gmsh",
        file,
        "-3",           # Generate a 3D mesh
        "-format",
        "msh2",
        "-o",
        mesh_filename
    ], capture_output=True, text=True)

    return mesh_filename   # Return the path to the generated mesh file


def parse_mesh(file: str) -> None:
    # This function parses a Gmsh mesh file to extract the nodes and elements (triangles and quadrilaterals). 
    # It reads the mesh file, identifies the sections for nodes and elements, and stores the relevant information 
    # in numpy arrays. Finally, it saves the parsed data to a .mat file for later use.
    with open(file, 'r') as f:
        lines = f.readlines()               # Read all lines from the mesh file into a list

        i = 0
        while i < len(lines):               # Loop through the lines using an index to keep track of the current position   
            line = lines[i].strip()
            if line == "$Nodes":            # Handle nodes
                n_nodes = int(lines[i + 1].strip())                 # Number of nodes
                nodes = np.zeros((n_nodes, 3), dtype=np.float64)    # Store all nodes
                for j in range(n_nodes):                            # Get all node coordinates
                    parts = lines[i + 2 + j].split()
                    x_coord = float(parts[-3])
                    y_coord = float(parts[-2])
                    z_coord = float(parts[-1])
                    nodes[j] = (x_coord, y_coord, z_coord)
                i += n_nodes + 1            # Move the index past the nodes section

            elif line == "$Elements":      # Handle elements
                n_elements = int(lines[i + 1].strip())      # Number of elements
                triangles = []
                outer_faces = []

                for j in range(n_elements):
                    parts = lines[i + 2 + j].split()
                    elem_type = int(parts[1])               # Get element type

                    # gmsh triangle = type 2
                    if elem_type == 2:                      # Check if the element is a triangle (type 2 in Gmsh)
                        triangles.append((np.float64(parts[-3]), 
                                          np.float64(parts[-2]), 
                                          np.float64(parts[-1])))
                        outer_faces.append((np.float64(parts[-3]), 
                                            np.float64(parts[-2]), 
                                            np.float64(parts[-1])))   # Add the triangle to the list of outer faces 
                    # gmsh quadrilateral = type 4
                    elif elem_type == 4:                    # Check if the element is a quadrilateral (type 4 in Gmsh)
                        p1 = np.float64(parts[-4])
                        p2 = np.float64(parts[-3])
                        p3 = np.float64(parts[-2])
                        p4 = np.float64(parts[-1])
                        faces = list(itertools.combinations([p1, p2, p3, p4], 3))   # Get all combinations of 3 nodes from the quadrilateral to form triangles
                        triangles.extend(faces)     # Add the triangles formed from the quadrilateral to the list of triangles
                triangles = np.array(triangles, dtype=np.float64)
                outer_faces = np.array(outer_faces, dtype=np.float64)
            i += 1

    savemat(f"{PARSED_DIR}{file[len(MESH_DIR):-4]}.mat", {            # Save the parsed nodes and triangles to a .mat file
        "nodes": nodes, 
        "triangles": triangles,
        "outer_faces": outer_faces   
        })

for (filename, refinement) in filenames:        # Loop through each geo file and refinement level, generate the mesh, and parse it
    generated_mesh = generate_meshes(filename, refinement)
    parse_mesh(generated_mesh)