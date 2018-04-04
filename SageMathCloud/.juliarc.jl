DIR=dirname(dirname(realpath(expanduser(("~/.juliarc.jl")))))
unshift!(LOAD_PATH, "$DIR/lib/julia")

# for Jupyter notebook on my screen
ENV["COLUMNS"] = 120

# I almost always want this, but don't error out if it
# break for e.g. https://github.com/timholy/Revise.jl/issues/77
try
    using Revise
end
