# Tsunami Simulation Lab Course

[Link](http://www5.in.tum.de/wiki/index.php/Bachelor-Praktikum:_Tsunami-Simulation_-_Summer_13)

## Dependencies

* [SCons](http://www.scons.org) (recommended for tests)
* [CxxTest](http://cxxtest.com) (for tests)
* [Doxygen](http://doxygen.org) (for documentation)

## Build Options

Build options can be supplied using 

    scons <command> [name=value [name=value] [...]]

Available options are

* `debug=1`: Use debug option `-g` (off by default)
* `optimize=<level>`: Use optimization level `-O<level>` (0 by default)
* `compiler=<name>`: Use different compiler, e.g. *clang* for LLVM/Clang

## Build

The default build command

    scons
    
will build the unit tests but not run it, useful for debugging propuses

## Tests

Tests can be run using

    scons check

## Documentation

HTML and LaTeX documentation can be generated using
    
    doxygen
    

