# About

This is a Julia CLI application to perform and visualize a 2D delaunay triangulation. It can be executed with Julia, or be compiled to a standalone application.

This is a learning project.


# Execution with Julia

To run the application with Julia, run `Triangulation.jl` in the Triangulation/src directory with Julia.

Print the arguments:

```
julia Triangulation.jl --help
```

Create the triangulation for a set of 100 random nodes:

```
julia Triangulation.jl --create-nodes=100
```


# Compilation of a standalone application

## Raw compilation

To compile into a standalone application, run `create_app.jl` in the Triangulation/build directory with Julia. This will create a directory named `exe_triangulation` next to the Triangulation directory. The binary and its libraries are in the bin sub-directory.

## Precompilation of the statements

To improve the startup time of the application, it is best to precompile the statements used. To do so:

1. Create a statements file by running `create_statements_file.jl` (in Triangulation/build) with julia. This will execute the functions of the project and list the statements used.
2. Run `build_app.jl` with Julia and add the `--precompile-statements` option. This will precompile the statements from the file generated at step 1.

## Further improvement on startup time

A strategy to reduce the startup time consists in precompiling as many of statements used by Julia as possible.

In the previous section, an attempt has been made to find all the statements behind the functions of the project. However, Julia may execute additional statements at startup.

Here is a strategy to trace and precompile those additional statements:

1. Start a Julia session from the system image created at the previous section. It is a file named `sys.dll` in the bin directory of a build. The compilation traces are activated and recorded to a file :
```
julia --trace-compile=add_statements_file.jl -J sys.dll
```

2. If `add_statements_file.jl` is not empty, append its content to the statements file created previously (`statements_file.jl` in the `build` subdirectory).

3. Delete or rename the executable, and perform again the build with precompilation of the statements (step 2 in the previous section).


# Motivation

This project is built in the process of learning Julia. It focuses on:
- the creation of a CLI application using Julia
- the compilation to a standalone application
- strategies to improve the startup time


# Version

This project was developed and tested with Julia 1.6.

