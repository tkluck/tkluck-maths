check:
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 test/julia/MFDeformations.jl
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 test/julia/Modules.jl
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 test/julia/Groebner.jl
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 test/julia/PolynomialRings.jl

mf:
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 ./matrices.jl

push:
	rsync -av --delete . 56861fef446d47f19213c4f066cb5ee4@compute4-us.sagemath.com:~/tkluck-maths

time:
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 ./time.jl
