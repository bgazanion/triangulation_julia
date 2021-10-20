# call this script from command line to compile
#   julia compile.jl [opt]
# optional: add --precompile to precompile the statements from statements_file.jl

using PackageCompiler
using Base.Filesystem

function build(args)
    cd("../..")

    # options
    # --precompile: precompile the statements using an existing statements file
    precompile_statements = false
    for arg_str in args
        if arg_str == "--precompile-statements"
            precompile_statements = true
        else
            println("Error: invalid option "*arg_str)
            exit()
        end
    end

    if precompile_statements
        create_app("Triangulation", 
            "exe_triangulation_precompiled",
            precompile_statements_file="./Triangulation/build/statements_file.jl")
    else
        create_app("Triangulation", 
        "exe_triangulation")
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    build(ARGS)
end