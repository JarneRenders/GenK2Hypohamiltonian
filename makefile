compiler=gcc
flags=-std=gnu11 -march=native
densenauty32=-DWORDSIZE=32 -DMAXN=WORDSIZE nautyW1.a
densenauty64=-DWORDSIZE=64 -DMAXN=WORDSIZE nautyL1.a

# The 32-bit version of this program takes less memory but only supports graphs up to 32 vertices.
32bit: genK2Hypohamiltonian.c bitset.h hamiltonicityMethods.c
	$(compiler) -o genK2Hypohamiltonian genK2Hypohamiltonian.c hamiltonicityMethods.c planarity.c $(densenauty32) $(flags) -O3

# The 64-bit version of this program.
64bit: genK2Hypohamiltonian.c bitset.h hamiltonicityMethods.c
	$(compiler) -o genK2Hypohamiltonian-64 genK2Hypohamiltonian.c hamiltonicityMethods.c planarity.c $(densenauty64) $(flags) -O3

profile: genK2Hypohamiltonian.c bitset.h hamiltonicityMethods.c
	$(compiler) -g -pg -o genK2Hypohamiltonian-pr genK2Hypohamiltonian.c hamiltonicityMethods.c planarity.c $(densenauty32) $(flags)

# Only execute until the common part of parallel parts is finished.
countsplit: genK2Hypohamiltonian.c hamiltonicityMethods.c bitset.h
	$(compiler) -DCOUNT_SPLIT_TIME -o genK2Hypohamiltonian-cs genK2Hypohamiltonian.c hamiltonicityMethods.c planarity.c $(densenauty32) $(flags) -O3

all: 32bit 64bit

.PHONY: clean
clean:
	rm -f genK2Hypohamiltonian genK2Hypohamiltonian-64 genK2Hypohamiltonian-cs genK2Hypohamiltonian-pr