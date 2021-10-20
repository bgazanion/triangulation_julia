# This file defines the elements required to create an app with PackageCompiler
# - module
# - julia_main function called at the execution of the app

module Triangulation

include("lib.jl")

# disable logging for levels Info and below
disable_logging(Logging.Info)

function julia_main()::Cint
    try
        main(ARGS)
        return 0
    catch
        return 1
    end
end

end # module

# set log level to WARNING for the module

#ENV["JULIA_WARNING"] = Triangulation

# CLI call from a terminal: call the main function with CLI args
if abspath(PROGRAM_FILE) == @__FILE__
    Triangulation.julia_main()
end
;