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
