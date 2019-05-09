# srep_init
Standalone srep initializer program based on libigl

Still under very active development (API and behavior subject to change)

## Install instructions
Clone this repo with

`git clone git@github.com:meglstuart/srep_init.git`

Download libigl with

`git clone --recurse-submodules https://github.com/libigl/libigl.git`

If you clone srep_init to `/path/to/srep_init` a good place to install libigl is `/path/to/libigl`

## Build
~~~
cd /path/to/srep_init
mkdir build
cd build
cmake ..
make
~~~

## Usage
Use `./srep_init --help` for full listing of options.
There is a (pretty incomplete) GUI that you can use with `./srep_init -g`, but for most usability I recommend the CLI

The generated srep is in `<output_folder>/initial_srep/` and the input mesh transformed to the same coordinate system is in `<output_folder>/aligned_mesh_for_refinement.vtk`.
