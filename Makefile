check:
	JULIA_LOAD_PATH="`pwd`/lib" julia -O3 test.jl

try:
	JULIA_LOAD_PATH="`pwd`/lib" julia -O3 mf-deformation.jl
