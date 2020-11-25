#CXX_FLAGS := -g

c2c_len3125: c2c_len3125.cpp
	/opt/rocm/bin/hipcc ${CXX_FLAGS} -I/workspace/rocfft/build/release/rocfft-install/include/ -o $@ $^
	extractkernel -i ./$@
	vim $@-000-gfx906.isa

