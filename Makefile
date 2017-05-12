check:
	JULIA_LOAD_PATH="`pwd`/lib" /home/tkluck/.local/bin/julia -O3 test/MFDeformations.jl
	JULIA_LOAD_PATH="`pwd`/lib" /home/tkluck/.local/bin/julia -O3 test/Modules.jl
	JULIA_LOAD_PATH="`pwd`/lib" /home/tkluck/.local/bin/julia -O3 test/Groebner.jl
	JULIA_LOAD_PATH="`pwd`/lib" /home/tkluck/.local/bin/julia -O3 test/PolynomialRings.jl

try:
	JULIA_LOAD_PATH="`pwd`/lib" /home/tkluck/.local/bin/julia -O3 ./matrices.jl
