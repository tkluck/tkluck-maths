check:
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 test/julia/GröbnerSingular.jl
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 test/julia/Coefficients.jl
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 test/julia/MFDeformations.jl
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 test/julia/Modules.jl

mf:
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 ./matrices.jl

push:
	git push bitbucket
	git push github
	rsync -av --delete . 56861fef446d47f19213c4f066cb5ee4@ssh.cocalc.com:~/tkluck-maths

time:
	JULIA_LOAD_PATH="`pwd`/lib/julia" julia -O3 ./time.jl
