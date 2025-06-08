# altel_acts

## Introduction


## Compiling


```sh
git clone https://github.com/eyiliu/altel_acts.git
cd altel_acts
mkdir BUILD
cd BUILD
cmake -DALTEL_EXTERNAL_CONFIG_ONLY=ON ..
make external_acts
make external_eudaq
make external_glfw
cmake ..
make
make install
cd ../INSTALL
```
The cmake will require net to download the acts, eudaq and glfw code from remote  git repo and checkout coressponding commit tags.

CMAKE_INSTALL_PREFIX is set to ** INSTALL **



## Using 


## System requirements
altel_acts is written and tesed on EL8 (RockyLinux8)


## Dependencies

+ eigen3-devel
+ boost169-devel
+ msgpack-devel
  
+ libX11-devel
+ libXpm-devel
+ libXft-devel
+ libXext-devel
+ mesa-libGL-devel
+ mesa-libGLU-devel
+ glew-devel
+ ftgl-devel
+ libXi-devel
+ libXi-devel
+ libXinerama-devel
+ glfw-devel 


## Changelog
