DIR=dirname(dirname(realpath(expanduser(("~/.juliarc.jl")))))
unshift!(LOAD_PATH, "$DIR/lib/julia")
