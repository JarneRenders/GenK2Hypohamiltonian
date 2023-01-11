# GenK2Hypohamiltonian
This repository contains a generator for K2-hypohamiltonian graphs which is created for the article "J. Goedgebeur, J. Renders and C.T. Zamfirescu, Generation and New Infinite Families of K<sub>2</sub>-hypohamiltonian Graphs, manuscript".

The latest version of this program can be obtained from <https://github.com/JarneRenders/GenK2Hypohamiltonian>.

This program can be used to generate all pairwise non-isomorphic K<sub>2</sub>-hypohamiltonian graphs of a specified order. It makes use of datastructures and methods from [`nauty`](https://pallini.di.uniroma1.it/) and [`K2-Hamiltonian Graphs`](https://github.com/JarneRenders/K2-Hamiltonian-Graphs).

### Installation

This requires a working shell and `make`.

- Download, extract and configure [`nauty`](https://pallini.di.uniroma1.it/).
- Compile the nauty libraries using: `make nautyW1.a` and `make nautyL1.a`.
- Copy the following files to the directory containing genK2Hypohamiltonian.c: 
	* naurng.h, 
	* nausparse.h, 
	* naututil.h, 
	* nauty.h, 
	* nautyL1.a, 
	* nautyW1.a, 
	* planarity.c, 
	* planarity.h, 
	* splay.c
- Compile using: 
	* `make` to create a binary for the 32 bit version
	* `make 64bit` to create a binary for the 64 bit version
	* `make all` to create both of the above
	* `make countsplit` to create a binary that stops execution as soon as the common part of each parallell part is done. Can be used to see the effect of the splitlevel on the runtime common part. 

The 32 bit version is significantly faster than the 64 bit version and allows the generation of Mandatory arguments to long options are mandatory for short options too.
 up to 32 vertices. In the rare cases when it is possible to generate K<sub>2</sub>-hypohamiltonian graphs on higher orders the 64 bit version allows the generation of K<sub>2</sub>-hypohamiltonian graphs up to 64 vertices. Use `make clean` to remove all binaries created in this way.

### Usage of GenK2Hypohamiltonian

This helptext can be found by executing `./genK2Hypohamiltonian -h`.

Usage: `./genK2Hypohamiltonian [-b] [-p] [-g#] [-d#] [-D#] [-X#] [-e] [-h] n [res/mod]`

Generate all pairwise non-isomorphic K<sub>2</sub>-hypohamiltonian graphs of order `n`. Graphs are sent to stdout in graph6 format. For more information on the format, see <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>.

The `res/mod` argument, should always appear after the specified order `n`. Otherwise, the order in which the arguments appear does not matter. Be careful not to put an argument immediately after one with an option. E.g. -g#b will not recognise the -b argument.

Mandatory arguments to long options are mandatory for short options too.
```
  -b, --bipartite           only generate bipartite K2-hypohamiltonian
                             graphs
  -p, --planar              only generate planar K2-hypohamiltonian graphs
  -g, --girth=GIRTH         only generate K2-hypohamiltonian graphs of
                             girth at least GIRTH
  -d, --minimum-degree=DEG  only generate K2-hypohamiltonian graphs with
                             minimum degree at least DEG
  -D, --maximum-degree=DEG  only generate K2-hypohamiltonian graphs with
                             maximum degree at most DEG
  -X, --splitlevel=LVL      set the splitlevel to LVL; used for
                             parallellization; a higher level means more
                             uniformity in running times across parts at 
                             the cost of longer runtime
  -e, --edges-count         print table of how many intermediate graphs had
                             a certain size
  -h, --help                print help message
  res/mod                   split the generation in mod (not necessarily
                             equally big) parts. Here part res will be 
                             executed. Splitting can cause a lot of extra
                             overhead and duplicate graphs. Make sure to
                             filter isomorphic graphs afterwards!
```

### Examples
`./genK2Hypohamiltonian 15`
Generates all pairwise non-isomorphic K<sub>2</sub>-hypohamiltonian graphs of order 15.

`./genK2Hypohamiltonian -p 15`
Generates all pairwise non-isomorphic planar K<sub>2</sub>-hypohamiltonian graphs of order 15.

`./genK2Hypohamiltonian -D3 15`
Generates all pairwise non-isomorphic K<sub>2</sub>-hypohamiltonian graphs of order 15 with maximum degree at most 3.

`./genK2Hypohamiltonian -d3 -D3 15`
Generates all pairwise non-isomorphic cubic K<sub>2</sub>-hypohamiltonian graphs of order 15.

`./genK2Hypohamiltonian -g5 -d3 -D3 15`
Generates all pairwise non-isomorphic cubic K<sub>2</sub>-hypohamiltonian graphs of order 15 and girth at least 5.

`./genK2Hypohamiltonian -p -g5 -d3 -D3 15`
Generates all pairwise non-isomorphic planar cubic K<sub>2</sub>-hypohamiltonian graphs of order 15 and girth at least 5.

`./genK2Hypohamiltonian 15 0/8`
Generates pairwise non-isomorphic K<sub>2</sub>-hypohamiltonian graphs of order 15 found in part 0. Doing this for parts 0-7 gives us all such graphs. The different parts can have isomorphic graph.