# About

This is a Julia CLI application to perform and visualize a 2D delaunay triangulation. It can be executed with Julia, or be compiled to a standalone application.


# Execution with Julia

To run the application with Julia, run Triangulation.jl in the Triangulation/src directory with Julia.

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

To compile into a standalone application, run create_app.jl in the Triangulation/build directory with Julia. This will create a directory named `exe_triangulation` next to the Triangulation directory. The binary and its libraries are in the bin sub-directory.

## Precompilation of the statements

To improve the startup time of the application, it is best to precompile the statements used. To do so:

1. Create a statements file by running create_statements_file.jl (in Triangulation/build) with julia. This will execute the functions of the project and list the statements used.
2. Run build_app.jl with Julia and add the --precompile-statements option. This will precompile the statements from the file generated at step 1.


# Motivation

This project is built in the process of learning Julia. It focuses on:
- the creation of a CLI application using Julia
- the compilation to a standalone application
- strategies to improve the startup time


