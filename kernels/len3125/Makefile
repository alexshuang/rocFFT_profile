#CXX_FLAGS := -g

c2c_len3125: c2c_len3125.cpp
	/opt/rocm/bin/hipcc ${CXX_FLAGS} -I/workspace/rocfft/build/release/rocfft-install/include/ -o $@ $^
	extractkernel -i ./$@
	../scripts/strip_isa.sh $@-000-gfx906.isa
	#vim $@-000-gfx906.isa -O rocfft_kernel_3125.h

