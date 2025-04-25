
# MeshSubdivider for Nuke

This Nuke plugin integrates OpenSubdiv to generate a smooth mesh from input geometry. OpenSubdiv is an open-source subdivision surface library developed by Pixar, designed for high-performance subdivision surface evaluation. 


![subdivider-AVLz247jO6C5kLvG](https://github.com/user-attachments/assets/4b224fce-67c0-4e6e-a9b2-2ee60c967eec)


## Features

Subdivision Scheme Control:

Scheme Type: This knob lets you choose the subdivision scheme.

Bilinear: Performs linear interpolation.

Catmull-Clark: Applies Catmull-Clark subdivision, which inserts a vertex at the face and edge midpoints and averages them. This is the scheme used in Pixar's RenderMan.

Boundary Interpolation:

This knob controls how boundaries (edges and vertices that are only connected to one face) are handled during subdivision. The options are:

boundary_none: No special boundary treatment.

boundary edge only: Only boundary edges are interpolated.

boundary edge and corner: Both boundary edges and corner vertices are interpolated.


## Installation


```bash
export NUKE_INCLUDE_DIR=/nuke_path/include
export NUKE_LIB_DIR=/nuke_path
mkdir build && cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=/vcpkg_path/scripts/buildsystems/vcpkg.cmake
make
```
    
## Known Issues

UV Interpolation: UVs are not perfectly interpolated by the OpenSubdiv node. This is a known issue and will be addressed in a future update. (TODO)
## License

This Nuke plugin is released under the [MIT](https://choosealicense.com/licenses/mit/) License. See LICENSE.txt for more information.
