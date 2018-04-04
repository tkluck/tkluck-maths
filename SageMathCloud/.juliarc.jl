DIR=dirname(dirname(realpath(expanduser(("~/.juliarc.jl")))))
unshift!(LOAD_PATH, "$DIR/lib/julia")

# for Jupyter notebook on my screen
ENV["COLUMNS"] = 120
