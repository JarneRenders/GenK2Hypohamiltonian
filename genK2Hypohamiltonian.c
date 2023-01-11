/**
 * genK2Hypohamiltonian.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 * 
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "bitset.h"
#include "hamiltonicityMethods.h"
#include "planarity.h" //Includes nauty.h

#define USAGE \
"\nUsage:./genK2Hypohamiltonian [-b] [-p] [-g#] [-d#] [-D#] [-X#] [-e] [-h] n [res/mod]\n\n"

#define HELPTEXT \
"Generate all pairwise non-isomorphic K2-hypohamiltonian graphs of order `n`.\n\
Graphs are sent to stdout in graph6 format. For more information on the format,\n\
see http://users.cecs.anu.edu.au/~bdm/data/formats.txt.\n\
\n\
The `res/mod` argument, should always appear after the specified order `n`.\n\
Otherwise, the order in which the arguments appear does not matter. Be careful\n\
not to put an argument immediately after one with an option. E.g. -g#b will not\n\
recognise the -b argument.\n\
\n\
Mandatory arguments to long options are mandatory for short options too.\n\
    -b, --bipartite             only generate bipartite K2-hypohamiltonian\n\
                                 graphs\n\
    -p, --planar                only generate planar K2-hypohamiltonian graphs\n\
    -g, --girth=GIRTH           only generate K2-hypohamiltonian graphs of\n\
                                 girth at least GIRTH\n\
    -d, --minimum-degree=DEG    only generate K2-hypohamiltonian graphs with\n\
                                 minimum degree at least DEG\n\
    -D, --maximum-degree=DEG    only generate K2-hypohamiltonian graphs with\n\
                                 maximum degree at most DEG\n\
    -X, --splitlevel=LVL        set the splitlevel to LVL; used for\n\
                                 parallellization; a higher level means more\n\
                                 uniformity in running times across parts at \n\
                                 the cost of longer runtime\n\
    -e, --edges-count           print table of how many intermediate graphs had\n\
                                 a certain size\n\
    -h, --help                  print help message\n\
    res/mod                     split the generation in mod (not necessarily\n\
                                 equally big) parts. Here part res will be \n\
                                 executed. Splitting can cause a lot of extra\n\
                                 overhead and duplicate graphs. Make sure to\n\
                                 filter isomorphic graphs afterwards!\n"

//  Used in parallellisation
#define SPLITLEVEL 1

//  Macro's for nauty representation
#define FOREACH(element,nautySet)\
 for(int element = nextelement((nautySet),MAXM,-1); (element) >= 0;\
 (element) = nextelement((nautySet),MAXM,(element)))
#define REMOVEONEEDGE(g,i,j,MAXM)\
 DELELEMENT(GRAPHROW(g,i,MAXM),j); DELELEMENT(GRAPHROW(g,j,MAXM),i)

// Graph structure containing nauty and bitset representations.
struct graph {
    graph nautyGraph[MAXN];
    int numberOfVertices;
    int numberOfEdges;
    bitset* adjacencyList;
    bitset* forbiddenEdges;
    bitset* verticesOfDeg;
}; 

//  Struct for passing options.
struct options {
    int minimalGirth;
    int minimumDegree;
    int maximumDegree;
    bool planarFlag;
    bool bipartiteFlag;
    int modulo;
    int remainder;
    int splitLevel;
    int splitCounter;
};

//  Struct keeping counters
struct counters {
    long long unsigned int *nOfNonIsoGraphsWithEdges;
    long long unsigned int nOfGraphsFound;
    long long unsigned int nOfTimesCheckedGirth;
    long long unsigned int nOfTimesHadForbiddenGirth;
    long long unsigned int nOfTimesCheckedHamiltonicity;
    long long unsigned int nOfTimesWasHamiltonian;
    long long unsigned int nOfTimesCheckedIsomorphism;
    long long unsigned int nOfTimesWasIsomorphic;
    long long unsigned int nOfTimesCheckedTypeA;
    long long unsigned int nOfTimesContainedTypeA;
    long long unsigned int nOfTimesCheckedTypeB;
    long long unsigned int nOfTimesContainedTypeB;
    long long unsigned int nOfTimesCheckedTypeC;
    long long unsigned int nOfTimesContainedTypeC;
    long long unsigned int nOfTimesCheckedDegreeObstruction;
    long long unsigned int nOfTimesContainedDegreeObstruction;
    long long unsigned int nOfTimesCheckedTriangleObstruction;
    long long unsigned int nOfTimesContainedTriangleObstruction;
    long long unsigned int nOfTimesCheckedSquareObstruction;
    long long unsigned int nOfTimesContainedSquareObstruction;
    long long unsigned int nOfTimesCheckedArrowObstruction;
    long long unsigned int nOfTimesContainedArrowObstruction;
    long long unsigned int nOfTimesCheckedStarObstruction;
    long long unsigned int nOfTimesContainedStarObstruction;
    long long unsigned int nOfTimesNoObstructionApplied;
    long long unsigned int nOfTimesObstructionFound;
    long long unsigned int nOfTimesTypeAObstructionChosen;
    long long unsigned int nOfTimesTypeBObstructionChosen;
    long long unsigned int nOfTimesTypeCObstructionChosen;
    long long unsigned int nOfTimesDegreeObstructionChosen;
    long long unsigned int nOfTimesTriangleObstructionChosen;
    long long unsigned int nOfTimes4CycleObstructionChosen;
    long long unsigned int nOfTimesGeneral4CycleObstructionChosen;
    long long unsigned int nOfTimesGeneralTriangleObstructionChosen;
    long long unsigned int nOfTimesCheckedMaximumDegree;
    long long unsigned int nOfTimesMaximumDegreeExceeded;
    long long unsigned int nOfTimesCheckedPlanarity;
    long long unsigned int nOfTimesWasNonPlanar;
};

/************************************************************************************
 * 
 * 
 *                      Macro's for dealing with graphs
 * 
 * 
 ************************************************************************************/

//  Initializer for empty graph.
#define emptyGraph(g) EMPTYGRAPH((g)->nautyGraph, (g)->numberOfVertices, MAXM);\
 (g)->verticesOfDeg[0] = compl(EMPTY,(g)->numberOfVertices);\
 (g)->numberOfEdges = 0;\
 for(int i = 0; i < (g)->numberOfVertices; i++) {(g)->adjacencyList[i] = EMPTY;\
 (g)->forbiddenEdges[i] = EMPTY;} 

//  Add one edge. Assume that every edge has already deg2 or more.
#define addEdge(g,i,j) {ADDONEEDGE((g)->nautyGraph, (i), (j), MAXM);\
 add((g)->adjacencyList[i], j); add((g)->adjacencyList[j],i);\
 (g)->numberOfEdges++;\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[i]) - 1], i);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[i])], i);\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[j]) - 1], j);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[j])], j);\
}

//  Remove one edge.
#define removeEdge(g,i,j) {REMOVEONEEDGE((g)->nautyGraph, (i), (j), MAXM);\
 removeElement((g)->adjacencyList[i], j); removeElement((g)->adjacencyList[j],i);\
 (g)->numberOfEdges--;\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[i]) + 1], i);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[i])], i);\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[j]) + 1], j);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[j])], j);\
}

//  Add edge to g->forbiddenEdges.
#define forbidEdge(g,i,j) {add((g)->forbiddenEdges[i],j);\
 add((g)->forbiddenEdges[j],i);}

//  Remove edge from g->forbiddenEdges.
#define permitEdge(g,i,j) {removeElement((g)->forbiddenEdges[i],j);\
 removeElement((g)->forbiddenEdges[j],i);}

#define areNeighbours(g,i,j) contains((g)->adjacencyList[(i)], j)


/************************************************************************************
 * 
 * 
 *                      Definitions for Nauty's splay tree
 * 
 * 
 ************************************************************************************/
typedef struct SPLAYNODE {
    graph* canonForm;
    struct SPLAYNODE *left, *right, *parent;
} SPLAYNODE;

SPLAYNODE *splayTreeArray[MAXN*(MAXN-1)/2] = {NULL};

#define SCAN_ARGS

#define ACTION(p) 

#define INSERT_ARGS , graph gCan[], int numberOfVertices,\
 bool *isPresent

int compareSplayNodeToGraph(SPLAYNODE* p, graph gCan[], int numberOfVertices);

#define COMPARE(p) compareSplayNodeToGraph(p, gCan, numberOfVertices);

#define PRESENT(p) {(*isPresent) = true;}

#define NOT_PRESENT(p) {p->canonForm = gCan; (*isPresent) = false;}

#define LOOKUP_ARGS , graph gCan[], int numberOfVertices

#include "splay.c"


/************************************************************************************
 * 
 * 
 *                          Functions for printing data 
 * 
 * 
 ************************************************************************************/

void printGraph(struct graph *g) {
    for(int i = 0; i < g->numberOfVertices; i++) {
        fprintf(stderr, "%d:", i);
        FOREACH(neighbour, GRAPHROW(g->nautyGraph, i, MAXM)) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr,"\n");   
}

void printNautyGraph(graph g[], int numberOfVertices) {
    for(int i = 0; i < numberOfVertices; i++) {
        fprintf(stderr, "%d:", i);
        FOREACH(neighbour, GRAPHROW(g, i, MAXM)) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr,"\n");   
}

void printAdjacencyList(struct graph *g) {
    for(int i = 0; i < g->numberOfVertices; i++) {
        fprintf(stderr, "%d:", i);
        forEach(neighbour, g->adjacencyList[i]) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

void printForbiddenEdges(struct graph *g) {
    for(int i = 0; i< g->numberOfVertices; i++) {
        forEach(nbr, g->forbiddenEdges[i]) {
            fprintf(stderr, " (%d,%d)\n", i, nbr);
        }
    }
}

//  Print gCan to stdout in graph6 format. 
void writeToG6(graph gCan[], int numberOfVertices) {
    char graphString[8 + numberOfVertices*(numberOfVertices - 1)/2];
    int pointer = 0;

    //  Save number of vertices in the first one, four or 8 bytes.
    if(numberOfVertices <= 62) {
        graphString[pointer++] = (char) numberOfVertices + 63;
    }
    else if(numberOfVertices <= 258047) {
        graphString[pointer++] = 63 + 63;
        for(int i = 2; i >= 0; i--) {
            graphString[pointer++] = (char) (numberOfVertices >> i*6) + 63;
        }
    }
    else if(numberOfVertices <= 68719476735) {
        graphString[pointer++] = 63 + 63;
        graphString[pointer++] = 63 + 63;
        for(int i = 5; i >= 0; i--) {
            graphString[pointer++] = (char) (numberOfVertices >> i*6) + 63;
        }
    }
    else {
        fprintf(stderr, "Error: number of vertices too large.\n");
        exit(1);
    }

    // Group upper triangle of adjacency matrix in groups of 6. See B. McKay's 
    // graph6 format.
    int counter = 0;
    char charToPrint = 0;
    for(int i = 1; i < numberOfVertices; i++) {
        for(int j = 0; j < i; j++) {
            charToPrint = charToPrint << 1;
            if(ISELEMENT(GRAPHROW(gCan, i, MAXM), j)) {
                charToPrint |= 1;
            }
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
                charToPrint = 0;
                counter = 0;
            }
        }
    }

    //  Pad final character with 0's.
    if(counter != 0) {
        while(counter < 6) {
            charToPrint = charToPrint << 1;
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
            }
        }
    }

    //  End with newline and end of string character.
    graphString[pointer++] = '\n';
    graphString[pointer++] = '\0';
    printf("%s", graphString);
}

/************************************************************************************/

//  Uses Nauty to get canonical form.
void createCanonicalForm(graph g[], graph gCan[], int numberOfVertices) {

    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    statsblk stats;

    densenauty(g, lab, ptn, orbits, &options, &stats, MAXM,
     numberOfVertices, gCan);
}

//  Splaynode contains the canonical form of a graph checked earlier.
int compareSplayNodeToGraph(SPLAYNODE* p, graph gCan[], int numberOfVertices) {
    return memcmp(p->canonForm, gCan, numberOfVertices * sizeof(graph));
}

//  Check recursively whether a path can be extended to a cycle of length
//  smaller than minimalGirth.
bool canBeForbiddenCycle(struct graph *g, bitset nbrsOfStart, int minimalGirth,
 int secondToLast, int last, int pathLength) {
    if(!isEmpty(intersection(nbrsOfStart, g->adjacencyList[last]))
     && pathLength + 1 < minimalGirth) {
        return true;
    }
    if(pathLength + 1 >= minimalGirth) {
        return false;
    }
    forEach(neighbourOfLast,
     difference(g->adjacencyList[last], singleton(secondToLast))) {
        if(canBeForbiddenCycle(g, nbrsOfStart, minimalGirth, last,
         neighbourOfLast, pathLength + 1)) {
            return true;
        }
    }
    return false;
}

// Check if there exists a cycle of length smaller than minimalGirth
// containing the specified edge.
bool containsForbiddenCycleWithEdge(struct graph *g, int minimalGirth,
 int endpoint1, int endpoint2) {
    bitset neighboursOfStartWithoutEndpoint2 =
     difference(g->adjacencyList[endpoint1], singleton(endpoint2));

    //  Loop over all neighbours of endpoint2 except endpoint1.
    forEach(neighbour,
     difference(g->adjacencyList[endpoint2], singleton(endpoint1))) {
        if(contains(g->adjacencyList[endpoint1], neighbour)
         && minimalGirth > 3) {
            return true;
        }

        //  Check if this path belongs to a forbidden cycle.
        if(canBeForbiddenCycle(g, neighboursOfStartWithoutEndpoint2,
         minimalGirth, endpoint2, neighbour, 3)) {
            return true;
        }
    }
    return false;
}

//  Method adapted from Nauty
int is_planar(struct graph *g) {

    if(g->numberOfEdges > (3*g->numberOfVertices - 6)) {
       return 0; //Can't be planar in this case
    }
    
    t_ver_sparse_rep V[g->numberOfVertices];
    t_adjl_sparse_rep A[2*g->numberOfEdges + 1];
    
    int k = 0;
    for(int i = 0; i < g->numberOfVertices; ++i)
        if(size(g->adjacencyList[i]) == 0)
            V[i].first_edge = NIL;
        else {
            V[i].first_edge = k;
            forEach(nbr, g->adjacencyList[i]) {
                A[k].end_vertex = nbr;
                A[k].next = k + 1;
                ++k;
            }
            A[k - 1].next = NIL;
        }

    if(k != 2 * g->numberOfEdges) {
        fprintf(stderr, "Error: decoding error (planarg)\n");
        exit(1);
    }
    
    int nbr_c;
    int edge_pos, v, w;
    t_dlcl **dfs_tree, **back_edges, **mult_edges;
    t_ver_edge *embed_graph;
    
    //TODO: possible optimisation: only start from last vertex?
    int ans = sparseg_adjl_is_planar(V, g->numberOfVertices, A, &nbr_c,
                                 &dfs_tree, &back_edges, &mult_edges,
                                 &embed_graph, &edge_pos, &v, &w);
    
    
    //Don't forget this otherwise out of mem!
    sparseg_dlcl_delete(dfs_tree, g->numberOfVertices);
    sparseg_dlcl_delete(back_edges, g->numberOfVertices);
    sparseg_dlcl_delete(mult_edges, g->numberOfVertices);
    embedg_VES_delete(embed_graph, g->numberOfVertices);
    
    return ans;
}

void extendWithEdge(struct graph *g, struct counters *counters, struct options *options, int
endpoint1, int endpoint2, int forbiddenEdges[], int *forbiddenEdgesLength);

//  Make edge forbidden and add to list forbiddenEdges, which will be used to
//  dynamically change the forbidden edges later.
void addForbiddenEdge(struct graph *g, int endpoint1, int endpoint2, int forbiddenEdges
[], int *forbiddenEdgesLength) {
    if(!contains(g->forbiddenEdges[endpoint1], endpoint2)) {
        forbiddenEdges[2*(*forbiddenEdgesLength)] = endpoint1;
        forbiddenEdges[2*(*forbiddenEdgesLength) + 1] = endpoint2;
        (*forbiddenEdgesLength)++;
        forbidEdge(g, endpoint1, endpoint2);
    }
}

//  Makes edges from list no longer forbidden. For dynamically keeping track
//  of forbidden edges.
void permitEdges(struct graph *g, int forbiddenEdges[], int forbiddenEdgesLength) {
    for(int i = 0; i < forbiddenEdgesLength; i++) {
        permitEdge(g,forbiddenEdges[2*i],forbiddenEdges[2*i+1]);
    }
}

//  Check whether adding this edge creates a bad cycle, hamiltonian graph or
//  isomorphic graph. If this is the case, we forbid the edge and add it to
//  the list forbiddenEdges. validEdges contains the edges that were checked
//  earlier by this function.
bool givesValidSuccessor(struct graph *g, struct counters *counters, struct options *options, int
forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], int endpoint1, int endpoint2) {
    if(contains(g->forbiddenEdges[endpoint1],endpoint2)) {
        return false;
    }
    if(contains(validEdges[endpoint1], endpoint2)) {
        return true;
    }

    addEdge(g, endpoint1, endpoint2);

    //  Check if this edge gives vertices a degree higher than the allowed maximum degree.
    if(options->maximumDegree != -1) {
        counters->nOfTimesCheckedMaximumDegree++;
        if(size(g->adjacencyList[endpoint1]) > options->maximumDegree ||
         size(g->adjacencyList[endpoint2]) > options->maximumDegree) {
                counters->nOfTimesMaximumDegreeExceeded++;
                if(!contains(g->forbiddenEdges[endpoint1], endpoint2)) {
                forbiddenEdges[2*(*forbiddenEdgesLength)] = endpoint1;
                forbiddenEdges[2*(*forbiddenEdgesLength) + 1] = endpoint2;
                (*forbiddenEdgesLength)++;
                forbidEdge(g, endpoint1, endpoint2);
            }
            removeEdge(g, endpoint1, endpoint2);
            return false;
        }
    }

    //  Check if girth is still at least minimalGirth.
    counters->nOfTimesCheckedGirth++;
    if(containsForbiddenCycleWithEdge(g, options->minimalGirth, endpoint1, endpoint2)) {
        counters->nOfTimesHadForbiddenGirth++;
        if(!contains(g->forbiddenEdges[endpoint1], endpoint2)) {
            forbiddenEdges[2*(*forbiddenEdgesLength)] = endpoint1;
            forbiddenEdges[2*(*forbiddenEdgesLength) + 1] = endpoint2;
            (*forbiddenEdgesLength)++;
            forbidEdge(g, endpoint1, endpoint2);
        }
        removeEdge(g, endpoint1, endpoint2);
        return false;
    }

    //  Check if hamiltonian.
    counters->nOfTimesCheckedHamiltonicity++;
    bitset remainingVertices = 
        compl(union(singleton(endpoint1), singleton(endpoint2)), g->numberOfVertices);
    if(canBeHamiltonian(g->adjacencyList, remainingVertices, endpoint2, endpoint1,
     g->numberOfVertices, 2)) {
        counters->nOfTimesWasHamiltonian++;
        if(!contains(g->forbiddenEdges[endpoint1], endpoint2)) {
            forbiddenEdges[2*(*forbiddenEdgesLength)] = endpoint1;
            forbiddenEdges[2*(*forbiddenEdgesLength) + 1] = endpoint2;
            (*forbiddenEdgesLength)++;
            forbidEdge(g, endpoint1, endpoint2);
        }
        removeEdge(g, endpoint1, endpoint2);
        return false;
    }

    removeEdge(g, endpoint1, endpoint2);
    add(validEdges[endpoint1], endpoint2);
    add(validEdges[endpoint2], endpoint1);
    return true;
}

//  Once the obstruction with the least number of successors is chosen, add
//  these successors (edges) recursively.
void addSuccessors(struct graph *g, struct counters *counters, struct options *options, int
newForbiddenEdges[], int *newForbiddenEdgesLength, int numberOfSuccessors, int *successors) {
    for(int i = 0; i < numberOfSuccessors; i++) {
        addEdge(g, successors[2*i], successors[2*i + 1]);
        extendWithEdge(g, counters, options, successors[2*i], successors[2*i + 1],
         newForbiddenEdges, newForbiddenEdgesLength);
        removeEdge(g, successors[2*i], successors[2*i + 1]);
    }
}

//  Assumes (W, V(G) - W) is a type A obs. For all edges working towards its
//  destruction, check how many are valid successors. Store corresponding
//  edges in successors[].
int countAllGoodAEdges(struct graph *g, struct counters *counters, struct options* options,
 int forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], bitset W, int successors[]) {
    int numberOfGoodAEdges = 0;

    // Add all good (W, V(G) - W) A-edges
    forEach(vertex, W) {
        forEachAfterIndex(newNbr,
         intersection(compl(g->adjacencyList[vertex], g->numberOfVertices), W), vertex) {
            if(givesValidSuccessor(g, counters, options, forbiddenEdges, forbiddenEdgesLength, 
             validEdges, vertex, newNbr)) {
                successors[2*numberOfGoodAEdges] = vertex;
                successors[2*numberOfGoodAEdges + 1] = newNbr;
                numberOfGoodAEdges++;
            }
        }
    }
    return numberOfGoodAEdges;
}

// int countAllGoodBEdges(struct graph *g, struct counters *counters, struct options* options, int
// forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], bitset setContainingTypeB, int
// successors[]) {
   
//     int numberOfGoodBEdges = 0;

//     //  Add all good B-edges.
//     forEach(vertex, setContainingTypeB) {
//         if(size(intersection(g->adjacencyList[vertex], setContainingTypeB)) <= 1) {
//             forEach(newNbr, intersection(compl(g->adjacencyList[vertex], g->numberOfVertices),
//              difference(setContainingTypeB, singleton(vertex)))) {
//                 if(givesValidSuccessor(g, counters, options, forbiddenEdges, forbiddenEdgesLength,
//                  validEdges, vertex, newNbr)) {
//                     successors[2*numberOfGoodBEdges] = vertex;
//                     successors[2*numberOfGoodBEdges + 1] = newNbr;
//                     numberOfGoodBEdges++;
//                 }
//             }
//         }
//     }
//     return numberOfGoodBEdges;
// }

// //  Compute k(G[W]).
// int computeKofW(struct graph *g, bitset W) {
//     if(g->numberOfVertices == 0) {
//         return 0;
//     }

//     int numberOfDeg1Vertices = 0;
//     int numberOfIsolatedVertices = 0;
//     int numberOfIsolatedK2s = 0;
//     bitset uncheckedVertices = W;
//     forEach(vertex, uncheckedVertices) {

//         //  Count isolated vertices and K2's.
//         int degreeInW = size(intersection(g->adjacencyList[vertex], W));
//         if(degreeInW == 0) {
//             numberOfIsolatedVertices++;
//         }
//         if(degreeInW == 1) {
//             numberOfDeg1Vertices++;
//             int neighbour = next(intersection(g->adjacencyList[vertex], W), -1);
//             if(size(intersection(g->adjacencyList[neighbour], W)) == 1) {
//                 numberOfDeg1Vertices++;
//                 numberOfIsolatedK2s++;
//             }
//             removeElement(uncheckedVertices, neighbour);
//         }
//     }
//     return numberOfIsolatedVertices + numberOfIsolatedK2s +
//              ((numberOfDeg1Vertices - 2*numberOfIsolatedK2s)/2 ? 1 : 0);
// }

//  Check whether (W,V(G)-W) is a type A obstruction under the assumption that
//  G[W] is a union of disjoint paths.
bool isTypeA(struct graph *g, bitset W, int numberOfDisjointPaths) {
    bitset X = compl(W, g->numberOfVertices);

    //  Check if X contains adjacent vertices.
    bool containsAdjacentVertices = false;
    forEach(vertex, X) {
        if(size(intersection(g->adjacencyList[vertex], X)) > 0) {
            containsAdjacentVertices = true;
            break;
        }
    }

    if(containsAdjacentVertices) {
        if(numberOfDisjointPaths >= size(X) - 1) {
            return true;
        }
    }

    return false;

}

//  Generate all vertex subsets W such that G[W] (induced subgraph) is a
//  union of disjoint paths of which none is a cycle.
void generateVertexDisjointPathSets(struct graph *g, struct counters *counters, struct
options *options, bitset W, bitset remainingVertices, int lastAdded, int numberOfDisjointPaths, int
forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], unsigned
int *numberOfLeastSuccessors, int leastSuccessors[]) {

    //  Check if current W is a type A obstruction.
    if(size(W) <= g->numberOfVertices - 3 && isTypeA(g, W, numberOfDisjointPaths)) {
        int successors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
        int numberOfSuccessors = countAllGoodAEdges(g, counters, options, forbiddenEdges,
         forbiddenEdgesLength, validEdges, W, successors);
        if(numberOfSuccessors < *numberOfLeastSuccessors) {
            *numberOfLeastSuccessors = numberOfSuccessors;
            memcpy(leastSuccessors, successors, 2*numberOfSuccessors*sizeof(int));
        }
    }

    //  The remainingVertices are only adjacent to vertices in the subgraph that
    //  have degree 0 or 1.
    forEachAfterIndex(vertex, remainingVertices, lastAdded) {
        int extraPaths = 0;

        //  Check if element to be added is adjacent to at most 1 element in the
        //  path. To generate all, this should be 2. Taking it to be 1 gives us
        //  no cycles which allows us to dynamically update it with practically
        //  no overhead.
        int degreeInPathCover = size(intersection(g->adjacencyList[vertex], W));
        if(degreeInPathCover > 1) {
            removeElement(remainingVertices, vertex);
            continue;
        }
        else if(degreeInPathCover == 0) {
            extraPaths++;
        }
        add(W, vertex);
        bitset remainingVerticesNew = remainingVertices;
        removeElement(remainingVerticesNew, vertex);

        //  The added vertex or its neighbours might have degree 2 in the
        //  W. If so, remove their neighbours as possible candidates to
        //  add to the subgraph. 
        if(size(intersection(g->adjacencyList[vertex], W)) > 1) {
            remainingVerticesNew = difference(remainingVerticesNew, g->adjacencyList[vertex]);
        }
        forEach(neighbour, intersection(g->adjacencyList[vertex], W)) {
            if(size(intersection(g->adjacencyList[neighbour], W)) > 1) {
                remainingVerticesNew = difference(remainingVerticesNew, g->adjacencyList[neighbour]);
            }
        }
        generateVertexDisjointPathSets(g, counters, options, W, remainingVerticesNew, vertex,
         numberOfDisjointPaths + extraPaths, forbiddenEdges, forbiddenEdgesLength, validEdges,
         numberOfLeastSuccessors, leastSuccessors);
        removeElement(W, vertex);
    }
}

//  leastSuccessors will contain all edges that belong to valid successors of
//  the type A obstruction which has the least valid successors.
bool getTypeAObstructionWithLeastSuccessors(struct graph *g, struct counters *counters, struct
options *options, int forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], unsigned
int *numberOfLeastSuccessors, int leastSuccessors[]) { 
    bitset W = EMPTY;
    bitset remainingVertices = compl(EMPTY, g->numberOfVertices);
    int numberOfDisjointPaths = 0;
    generateVertexDisjointPathSets(g, counters, options, W, remainingVertices, -1,
     numberOfDisjointPaths, forbiddenEdges, forbiddenEdgesLength, validEdges,
     numberOfLeastSuccessors, leastSuccessors);

    return *numberOfLeastSuccessors != -1;
}

//  Assumes (W, V(G) - W, vw) is a type C obs. For all edges working towards
//  its destruction, check how many are valid successors. Store corresponding
//  edges in successors[].
int countAllGoodCEdges(struct graph *g, struct counters *counters, struct options *options,
 int forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], bitset W,
 bitset goodCEdgeEndpointsInX, int successors[]) {
    int numberOfGoodCEdges = 0;
    forEach(vertexInW, W) {
        bitset nonNeighbours = difference(
            compl(g->adjacencyList[vertexInW], g->numberOfVertices),
            singleton(vertexInW));
        forEach(goodCNbr, intersection(nonNeighbours, union(W, goodCEdgeEndpointsInX))) {
            if(givesValidSuccessor(g, counters, options, forbiddenEdges, forbiddenEdgesLength,
             validEdges, vertexInW, goodCNbr)) {
                successors[2*numberOfGoodCEdges] = vertexInW;
                successors[2*numberOfGoodCEdges + 1] = goodCNbr;
                numberOfGoodCEdges++;
            }
        }
    }
    return numberOfGoodCEdges;
}

//  Check that there are v,w in V(G) - W such that (W, V(G)-W, vw) is a type C
//  obstruction under the assumption that W is an independent set of G and
//  store them in v and w respectively.
bool isTypeC(struct graph *g, bitset W, bitset *goodCEdgeEndpointsInX) {
    if(size(W) > 1) {

        //  Compute n1 and n2
        int n1 = 0;
        int n2 = 0;

        //  Used later to add good C edges
        bitset verticesInXWithOneOrLessNeighbourInW = EMPTY;
        bitset X = compl(W, g->numberOfVertices);

        //  Compute n1 and n2. Use this to also check whether G[X] contains
        //  adjacent vertices. We assume mindeg of G is 2.
        bool existsAdjacentPairInX = false;
        forEach(v, X) {
            if(size(intersection(g->adjacencyList[v], W)) == 1) {
                n1++;
                add(verticesInXWithOneOrLessNeighbourInW, v);
            }
            else if(size(intersection(g->adjacencyList[v], W)) > 1) {
                n2++;
                if(size(intersection(g->adjacencyList[v], X)) > 0) {
                    existsAdjacentPairInX = true;
                }
            }
            else {
                add(verticesInXWithOneOrLessNeighbourInW, v);
            }
        }
        //  Since mindeg is assumed to be 2, vertex with at most 1 neighbour in W has neighbour in X.
        if(size(verticesInXWithOneOrLessNeighbourInW) == 0 && !existsAdjacentPairInX) {
            return false;   
        }

        //  If this is the case there are no type C obstructions for this
        //  independent set.
        if(2*n2 + n1 < 2*size(W) + 4) {

            //  Find v,w s.t. (W,X,vw) is type C
            forEach(v, X) {
                int numberToSubtract = 0;
                if(size(intersection(g->adjacencyList[v], W)) == 1) {
                    numberToSubtract = 1; // v contributed to n1
                }
                else if(size(intersection(g->adjacencyList[v], W)) > 1) {
                    numberToSubtract = 2; // v contributed to n2
                }
                forEach(w, intersection(g->adjacencyList[v], W)) {
                    int numberToSubtractCopy = numberToSubtract;
                    if(size(intersection(g->adjacencyList[w], W)) == 1) {
                        numberToSubtract += 1; // w contributed to n1
                    }
                    else if(size(intersection(g->adjacencyList[w], W)) > 1) {
                        numberToSubtract += 2; // w contributed to n2
                    }
                    if(2*n2 + n1 - numberToSubtract < 2*size(W)) {

                        //  We store the vertices in X which have one or fewer
                        //  neighbours in W which are not v or w.
                        (*goodCEdgeEndpointsInX) = difference(verticesInXWithOneOrLessNeighbourInW,
                         union(singleton(v), singleton(w)));
                        return true;
                    }
                    numberToSubtract = numberToSubtractCopy;
                }
            }
        }
    }
    return false;
}

//  Generate all independent sets W using simple Bron Kerbosch algorithm.
//  Check if for some v,w in V(G) (W,W - V(G), vw) is a type C obs.
void findAllTypeCIndependentSets(struct graph *g, struct counters *counters, struct
options *options, bitset independentSet, bitset remainingVertices, int forbiddenEdges
[], int *forbiddenEdgesLength, bitset validEdges[], unsigned int *numberOfLeastSuccessors, int
leastSuccessors[]) {

    bitset goodCEdgeEndpointsInX;
    if(isTypeC(g, independentSet, &goodCEdgeEndpointsInX)) {
        int successors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
        int numberOfSuccessors = countAllGoodCEdges(g, counters, options, forbiddenEdges,
         forbiddenEdgesLength, validEdges, independentSet, goodCEdgeEndpointsInX, successors);
        if(numberOfSuccessors < *numberOfLeastSuccessors) {
            *numberOfLeastSuccessors = numberOfSuccessors;
            memcpy(leastSuccessors, successors, 2*numberOfSuccessors*sizeof(int));
        }
    }

    //BronKerbosch (simple)
    forEach(v, remainingVertices) {
        bitset independentSetNew = union(independentSet, singleton(v));
        bitset nonNeighboursOfV = compl(g->adjacencyList[v], g->numberOfVertices);
        removeElement(nonNeighboursOfV, v);
        bitset remainingVerticesNew = intersection(remainingVertices, nonNeighboursOfV);
        findAllTypeCIndependentSets(g, counters, options, independentSetNew, remainingVerticesNew,
         forbiddenEdges, forbiddenEdgesLength, validEdges, numberOfLeastSuccessors, leastSuccessors);
        removeElement(remainingVertices, v);
    }
}

//  leastSuccessors will contain all edges that belong to valid successors of
//  the type C obstruction which has the least valid successors. While I
//  believe number returned by this method will often be the least (I have
//  not checked this), in theory this is not guaranteed, since for a pair
//  (W,X), we stop as soon as we find v,w in X making (W,X,vw) a type C obs.
bool getTypeCIndepSetWithLeastSuccessors(struct graph *g, struct counters *counters, struct
options *options, int forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], unsigned
int *numberOfLeastSuccessors, int leastSuccessors[]) {  
    bitset independentSet = EMPTY;
    bitset remainingVertices = compl(EMPTY, g->numberOfVertices);
    findAllTypeCIndependentSets(g, counters, options, independentSet, remainingVertices,
     forbiddenEdges, forbiddenEdgesLength, validEdges, numberOfLeastSuccessors, leastSuccessors);

    return (*numberOfLeastSuccessors) != -1;
}

//  For many obstructions, the edges corresponding to the successors are
//  obtained by adding all v-edges with v some vertex except for maybe some
//  blockedVertices. We count how many of these are valid successors.
int countValidSuccessorsWithEndpoint(struct graph *g, struct counters *counters, struct
options *options, int endpoint, int successors[], bitset blockedVertices, int forbiddenEdges
[], int *forbiddenEdgesLength, bitset validEdges[]) {
    int numberOfGoodEdges = 0;
    bitset nonNeighbours = compl(g->adjacencyList[endpoint], g->numberOfVertices);
    removeElement(nonNeighbours, endpoint);
    forEach(v, nonNeighbours) {
        if(contains(blockedVertices, v)) {
            continue;
        }
        if(givesValidSuccessor(g, counters, options, forbiddenEdges, forbiddenEdgesLength,
         validEdges, endpoint, v)) {
            successors[2*numberOfGoodEdges] = endpoint;
            successors[2*numberOfGoodEdges + 1] = v;
            numberOfGoodEdges++;
        }
    }
    return numberOfGoodEdges;
}

//  Minimum degree can specified or should be a least three. 
bool getBadDegreeVertexWithLeastSuccessors(struct graph *g, struct counters *counters, struct
options *options, int forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], unsigned
int *numberOfLeastSuccessors, int leastSuccessors[]) {
    forEach(deg2Vertex, g->verticesOfDeg[2]) {
        int successors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
        int numberOfSuccessors = countValidSuccessorsWithEndpoint(g, counters, options, deg2Vertex,
         successors, EMPTY, forbiddenEdges, forbiddenEdgesLength, validEdges);
        if(numberOfSuccessors < *numberOfLeastSuccessors) {
            *numberOfLeastSuccessors = numberOfSuccessors;
            memcpy(leastSuccessors, successors, 2*numberOfSuccessors*sizeof(int));
        }
    }

    //  In case a minimum degree was chosen:
    for(int i = 3; i < options->minimumDegree; i++) {
        forEach(vertex, g->verticesOfDeg[i]) {
            int successors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
            int numberOfSuccessors = countValidSuccessorsWithEndpoint(g, counters, options, vertex,
             successors, EMPTY, forbiddenEdges, forbiddenEdgesLength, validEdges);
            if(numberOfSuccessors < *numberOfLeastSuccessors) {
                *numberOfLeastSuccessors = numberOfSuccessors;
                memcpy(leastSuccessors, successors, 2*numberOfSuccessors*sizeof(int));
            }
        }
    }
    return *numberOfLeastSuccessors != -1;
}

//  Any vertex on a triangle should have degree at least 4.
bool getDeg3InTriangleWithLeastSuccessors(struct graph *g, struct counters *counters, struct
options *options, int forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], unsigned
int *numberOfLeastSuccessors, int leastSuccessors[]) {
    forEach(vertex, g->verticesOfDeg[3]) {
        int previousNeighbour = next(g->adjacencyList[vertex], -1);
        forEachAfterIndex(nbr, g->adjacencyList[vertex], previousNeighbour) {
            if(!isEmpty(intersection(g->adjacencyList[vertex], g->adjacencyList[nbr]))) {
                int successors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
                int numberOfSuccessors = countValidSuccessorsWithEndpoint(g, counters, options,
                 vertex, successors, EMPTY, forbiddenEdges, forbiddenEdgesLength, validEdges);
                if(numberOfSuccessors < *numberOfLeastSuccessors) {
                    *numberOfLeastSuccessors = numberOfSuccessors;
                    memcpy(leastSuccessors, successors, 2*numberOfSuccessors*sizeof(int));
                }
            }
        }
    }
    return *numberOfLeastSuccessors != -1;
}

//  Any vertex on a 4-cycle should have degree at least 4.
bool getDeg3In4CycleWithLeastSuccessors(struct graph *g, struct counters *counters, struct
options *options, int forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], unsigned
int *numberOfLeastSuccessors, int leastSuccessors[]) {
    
    forEach(vertex, g->verticesOfDeg[3]) {
        forEach(firstNbr, g->adjacencyList[vertex]) {
            forEachAfterIndex(secondNbr, g->adjacencyList[vertex], firstNbr) {
                if(!isEmpty(intersection(difference(g->adjacencyList[firstNbr], singleton(vertex)),
                 difference(g->adjacencyList[secondNbr], singleton(vertex))))) {
                    int successors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
                    int numberOfSuccessors = countValidSuccessorsWithEndpoint(g, counters, options,
                     vertex, successors, EMPTY, forbiddenEdges, forbiddenEdgesLength, validEdges);
                    if(numberOfSuccessors < *numberOfLeastSuccessors) {
                        *numberOfLeastSuccessors = numberOfSuccessors;
                        memcpy(leastSuccessors, successors, 2*numberOfSuccessors*sizeof(int));
                    }
                }
            }
        }
    }
    return *numberOfLeastSuccessors != -1;
}

//  Look for (v,uvw)-obstructions, i.e. a vertex v on a triangle uvw with at
//  most one neighbour not in N(u) cup N(w).
bool getGeneralTriangleVertexWithLeastSuccessors(struct graph *g, struct counters *counters, struct
options *options, int forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], unsigned
int *numberOfLeastSuccessors, int leastSuccessors[]) {
    
    for(int i = 4; i < g->numberOfVertices - 1; i++) {
        forEach(v, g->verticesOfDeg[i]) {
            bool foundObstruction = false;
            forEach(u, g->adjacencyList[v]) {
                forEachAfterIndex(w, intersection(g->adjacencyList[v], g->adjacencyList[u]), u) {
                    if(areNeighbours(g, u, w)) { // uvw is a triangle.
                        if(size(intersection(g->adjacencyList[v],
                         union(g->adjacencyList[u], g->adjacencyList[w]))) >= i - 1) {
                            int successors[g->numberOfVertices*(g->numberOfVertices - 1) - 
                                2*g->numberOfEdges];
                            int numberOfSuccessors = countValidSuccessorsWithEndpoint(g, counters, 
                             options, v, successors, union(g->adjacencyList[u], g->adjacencyList[w]),
                             forbiddenEdges, forbiddenEdgesLength, validEdges);
                            if(numberOfSuccessors < *numberOfLeastSuccessors) {
                                *numberOfLeastSuccessors = numberOfSuccessors;
                                memcpy(leastSuccessors, successors, 2*numberOfSuccessors*sizeof(int));
                                foundObstruction = true;
                                break;
                            }
                        }
                    }
                    if(foundObstruction) {
                        break;
                    }
                }
                if(foundObstruction) {
                    break;
                }
            }
        }
    }
    return *numberOfLeastSuccessors != -1;
}

//  Look for (v,uvwx)-obstructions, i.e. a vertex v on a 4-cycle uvwx with at
//  most one neighbour not in N[x].
bool getGeneral4CycleVertexWithLeastSuccessors(struct graph *g, struct counters *counters, struct
options *options, int forbiddenEdges[], int *forbiddenEdgesLength, bitset validEdges[], unsigned
int *numberOfLeastSuccessors, int leastSuccessors[]) {
    
    for(int i = 4; i < g->numberOfVertices - 1; i++) {
        forEach(v, g->verticesOfDeg[i]) {
            bool foundObstruction = false;
            forEach(u, g->adjacencyList[v]) {
                forEachAfterIndex(w, g->adjacencyList[v], u) {
                    forEach(x, difference(intersection(g->adjacencyList[u], g->adjacencyList[w]),
                     singleton(v))) {
                        if(!areNeighbours(g, x, v)) {
                            if(size(intersection(g->adjacencyList[v], g->adjacencyList[x]))
                             >= i - 1) {
                                int successors[g->numberOfVertices*(g->numberOfVertices - 1)
                                 - 2*g->numberOfEdges];
                                int numberOfSuccessors = countValidSuccessorsWithEndpoint(g,
                                 counters, options, v, successors, g->adjacencyList[x],
                                 forbiddenEdges, forbiddenEdgesLength, validEdges);
                                if(numberOfSuccessors < *numberOfLeastSuccessors) {
                                    *numberOfLeastSuccessors = numberOfSuccessors;
                                    memcpy(leastSuccessors, successors,
                                     2*numberOfSuccessors*sizeof(int));
                                    foundObstruction = true;
                                    break;
                                }
                            }
                        }
                    }
                    if(foundObstruction) {
                        break;
                    }
                }
                if(foundObstruction) {
                    break;
                }
            }
        }
    }
    return *numberOfLeastSuccessors != -1;
}

void extendWithEdge(struct graph *g, struct counters *counters, struct options *options, int
endpoint1, int endpoint2, int forbiddenEdges[], int *forbiddenEdgesLength) {

    counters->nOfNonIsoGraphsWithEdges[g->numberOfEdges]++;

    //  This list will keep track of edges that became forbidden after adding
    //  the previous edge. Edges already forbidden should not be in the list.
    //  At the end of this iteration we will permit all edge of this list
    //  again.  We assume that our starting g has at least
    //  g->numberOfVertices + 3 edges.
    int newForbiddenEdges[(g->numberOfVertices)*(g->numberOfVertices-1)-2*(g->numberOfVertices + 3)];
    int newForbiddenEdgesLength = 0;

    //  Keeps track of non-edges deemed valid by isValidSuccessor. This way we
    //  do not do any double checking of edges (at the cost of some memory).
    bitset validEdges[g->numberOfVertices];
    for(int i = 0; i < g->numberOfVertices; i++) {
        validEdges[i] = EMPTY;
    }

    //  Due to the overhead of this method, we cannot check this too often.
    //  Hence, we check it here instead of in givesValidSuccessor(). Checking this
    //  every two times is much slower. 
    if(options->planarFlag) {
        counters->nOfTimesCheckedPlanarity++;
        if(!is_planar(g)) {
            counters->nOfTimesWasNonPlanar++;
            if(!contains(g->forbiddenEdges[endpoint1], endpoint2)) {
                forbiddenEdges[2*(*forbiddenEdgesLength)] = endpoint1;
                forbiddenEdges[2*(*forbiddenEdgesLength) + 1] = endpoint2;
                (*forbiddenEdgesLength)++;
                forbidEdge(g, endpoint1, endpoint2);
            }
            return;
        }
    }

    //  Check for isomorphism with an earlier intermediate graph.
    graph* gCan = malloc(sizeof(graph)*g->numberOfVertices);
    if(gCan == NULL) {
        fprintf(stderr, "Error: out of memory\n");
        exit(1);
    }
    createCanonicalForm(g->nautyGraph, gCan, g->numberOfVertices);
    bool isPresent;
    counters->nOfTimesCheckedIsomorphism++;
    splay_insert(&splayTreeArray[g->numberOfEdges], gCan, g->numberOfVertices, &isPresent);
    if(isPresent) {
        counters->nOfTimesWasIsomorphic++;
        if(!contains(g->forbiddenEdges[endpoint1], endpoint2)) {
            forbiddenEdges[2*(*forbiddenEdgesLength)] = endpoint1;
            forbiddenEdges[2*(*forbiddenEdgesLength) + 1] = endpoint2;
            (*forbiddenEdgesLength)++;
            forbidEdge(g, endpoint1, endpoint2);
        }
        free(gCan);
        return;
    }

    //  Used for paralellization. Checking the split needs to be done after
    //  the isomorphism check. The above check may forbid edges in the
    //  splitlevel-1 graphs. If we deal with the splitting prior to these
    //  checks, different remainders would give different forbidden edges in
    //  these graphs and hence the splitting counter would change depending
    //  on the choice of remainder.
    if(g->numberOfEdges - (g->numberOfVertices + 3) == options->splitLevel) {

        //  Count how long the common part of each split takes. Algorithm
        //  stops as soon as the common part is finished.
        #ifdef COUNT_SPLIT_TIME 
            return;
        #endif
        options->splitCounter++;
        if(options->splitCounter == options->modulo) {
            options->splitCounter = 0;
        }
        if(options->splitCounter != options->remainder) {
            permitEdges(g, newForbiddenEdges, newForbiddenEdgesLength);
            return;
        }
    }

    //  Find the obstruction with least number of valid successors. 
    unsigned int numberOfLeastSuccessors = -1;
    int *successors = NULL;
    char chosenObstruction;

    unsigned int numberOfTypeASuccessors = -1;
    int typeASuccessors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
    counters->nOfTimesCheckedTypeA++;
    bool hasTypeASuccessors = getTypeAObstructionWithLeastSuccessors(g, counters, options,
     newForbiddenEdges, &newForbiddenEdgesLength, validEdges, &numberOfTypeASuccessors, 
     typeASuccessors);
    if(hasTypeASuccessors) {
        counters->nOfTimesContainedTypeA++;
        if(numberOfTypeASuccessors < numberOfLeastSuccessors) {
            numberOfLeastSuccessors = numberOfTypeASuccessors;
            successors = typeASuccessors;
            chosenObstruction = 'A';
        }
    }

    unsigned int numberOfDegSuccessors = -1;
    int degSuccessors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
    bool hasBadDegVertex = getBadDegreeVertexWithLeastSuccessors(g, counters, options,
     newForbiddenEdges, &newForbiddenEdgesLength, validEdges, &numberOfDegSuccessors, 
     degSuccessors);
    counters->nOfTimesCheckedDegreeObstruction++;
    if(hasBadDegVertex) {
        counters->nOfTimesContainedDegreeObstruction++;
        if(numberOfDegSuccessors < numberOfLeastSuccessors) {
            numberOfLeastSuccessors = numberOfDegSuccessors;
            successors = degSuccessors;
            chosenObstruction = '2';
        }
    }

    unsigned int numberOfTriangleSuccessors = -1;
    int triangleSuccessors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
    if(options->minimalGirth < 4) {
        bool hasDeg3VertexInTriangle = getDeg3InTriangleWithLeastSuccessors(g, counters, options,
         newForbiddenEdges, &newForbiddenEdgesLength, validEdges, &numberOfTriangleSuccessors, 
         triangleSuccessors);
        counters->nOfTimesCheckedTriangleObstruction++;
        if(hasDeg3VertexInTriangle) {
            counters->nOfTimesContainedTriangleObstruction++;
            if( numberOfTriangleSuccessors < numberOfLeastSuccessors) {
                numberOfLeastSuccessors = numberOfTriangleSuccessors;
                successors = triangleSuccessors;
                chosenObstruction = '3';
            }
        }
    }

    unsigned int numberOf4CycleSuccessors = -1;
    int fourCycleSuccessors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
    if(options->minimalGirth < 5) {
        bool has4CycleObstruction = getDeg3In4CycleWithLeastSuccessors(g, counters, options,
         newForbiddenEdges, &newForbiddenEdgesLength, validEdges, &numberOf4CycleSuccessors, 
         fourCycleSuccessors);
        counters->nOfTimesCheckedSquareObstruction++;
        if(has4CycleObstruction) {
            counters->nOfTimesContainedSquareObstruction++;
            if( numberOf4CycleSuccessors < numberOfLeastSuccessors) {
                numberOfLeastSuccessors = numberOf4CycleSuccessors;
                successors = fourCycleSuccessors;
                chosenObstruction = '4';
            }
        }
    }

    unsigned int numberOfTypeCSuccessors = -1;
    int typeCSuccessors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
    bool hasTypeCObstruction = getTypeCIndepSetWithLeastSuccessors(g, counters, options, 
     newForbiddenEdges, &newForbiddenEdgesLength, validEdges, &numberOfTypeCSuccessors, 
     typeCSuccessors);
    counters->nOfTimesCheckedTypeC++;
    if(hasTypeCObstruction) {
        counters->nOfTimesContainedTypeC++;
        if( numberOfTypeCSuccessors < numberOfLeastSuccessors) {
            numberOfLeastSuccessors = numberOfTypeCSuccessors;
            successors = typeCSuccessors;
            chosenObstruction = 'C';
        }
    }

    unsigned int numberOfGeneralTriangleSuccessors = -1;
    int generalTriangle[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
    if(options->minimalGirth < 4) {
        bool hasGeneralTriangleObstruction = getGeneralTriangleVertexWithLeastSuccessors(g,
         counters, options, newForbiddenEdges, &newForbiddenEdgesLength, validEdges, 
         &numberOfGeneralTriangleSuccessors, generalTriangle);
        counters->nOfTimesCheckedStarObstruction++;
        if(hasGeneralTriangleObstruction) {
            counters->nOfTimesContainedStarObstruction++;
            if(numberOfGeneralTriangleSuccessors < numberOfLeastSuccessors) {
                numberOfLeastSuccessors = numberOfGeneralTriangleSuccessors;
                successors = generalTriangle;
                chosenObstruction = 's';
            }
        }
    }

    unsigned int numberOfGeneral4CycleSuccessors = -1;
    int general4CycleSuccessors[g->numberOfVertices*(g->numberOfVertices - 1) - 2*g->numberOfEdges];
    if(options->minimalGirth < 5) {
        bool hasGeneral4CycleObstruction = getGeneral4CycleVertexWithLeastSuccessors(g, counters,
         options, newForbiddenEdges, &newForbiddenEdgesLength, validEdges, 
         &numberOfGeneral4CycleSuccessors, general4CycleSuccessors);
        counters->nOfTimesCheckedArrowObstruction++;
        if(hasGeneral4CycleObstruction) {
            counters->nOfTimesContainedArrowObstruction++;
            if(numberOfGeneral4CycleSuccessors < numberOfLeastSuccessors) {
                numberOfLeastSuccessors = numberOfGeneral4CycleSuccessors;
                successors = general4CycleSuccessors;
                chosenObstruction = 'a';
            }
        }
    }

    switch(chosenObstruction) {
        case 'A':
            counters->nOfTimesTypeAObstructionChosen++;
            break;
        case 'C':
            counters->nOfTimesTypeCObstructionChosen++;
            break;
        case '2':
            counters->nOfTimesDegreeObstructionChosen++;
            break;
        case '3':
            counters->nOfTimesTriangleObstructionChosen++;
            break;
        case '4':
            counters->nOfTimes4CycleObstructionChosen++;
            break;
        case 's':
            counters->nOfTimesGeneralTriangleObstructionChosen++;
            break;
        case 'a':
            counters->nOfTimesGeneral4CycleObstructionChosen++;
            break;
    }

    //  Add all edges corresponding to valid successors of chosen
    //  obstruction.
    if(successors) {
        counters->nOfTimesObstructionFound++;
        addSuccessors(g, counters, options, newForbiddenEdges, &newForbiddenEdgesLength,
         numberOfLeastSuccessors, successors);
        permitEdges(g, newForbiddenEdges, newForbiddenEdgesLength);
        return;
    }

    //  We only get here if no obstructions are found.
    counters->nOfTimesNoObstructionApplied++;
    int temp[2] = {-1,-1};
    if(isK2Hamiltonian(g->adjacencyList, g->numberOfVertices, false, false,
     temp)) {
        writeToG6(gCan, g->numberOfVertices);
        counters->nOfGraphsFound++;
    }

    //  If no obs found, add all valid edges.
    for(int i = 0; i < g->numberOfVertices - 2; i++) {
        bitset nonNeighbours = compl(g->adjacencyList[i], g->numberOfVertices);
        removeElement(nonNeighbours, i);
        forEach(j, nonNeighbours) {
            if(contains(g->forbiddenEdges[i],j)) {
                continue;
            }
            if(givesValidSuccessor(g, counters, options, forbiddenEdges, forbiddenEdgesLength,
             validEdges, i, j)) {
                addEdge(g, i, j);
                extendWithEdge(g, counters, options, i, j, newForbiddenEdges, &newForbiddenEdgesLength);
                removeEdge(g, i, j);
            }
        }
    }
    permitEdges(g, newForbiddenEdges, newForbiddenEdgesLength);
}


//  Initialization of the graph with a cycle and a K2. We also add the first
//  four edges.
void addK2Edges(struct graph *g, struct counters *counters, struct options
*options) {

    //  First edge always connects to 0 up to isomorphism.
    addEdge(g, g->numberOfVertices - 2, 0);

    //  Second edge goes from n-1 to the cycle. It cannot join the cycle next to
    //  0, since that would give a hamiltonian cycle. Depending on the girth
    //  restriction it cannot join the cycle girth - 4 places next to 0.
    bitset possibilitiesForSecondEdge = 
     compl(EMPTY, g->numberOfVertices - 2);
    removeElement(possibilitiesForSecondEdge, g->numberOfVertices - 3);
    removeElement(possibilitiesForSecondEdge, 1);

    //  We also forbid these edges since no supergraph should contain them.
    forbidEdge(g, g->numberOfVertices - 1, g->numberOfVertices - 3);
    forbidEdge(g, g->numberOfVertices - 1, 1);

    //  If maximum degree was chosen, we cannot allow vertices of degree D+1.
    //  However, we can at most have vertices with degree 4 in this method.
    if(options->maximumDegree == 3) {
        removeElement(possibilitiesForSecondEdge, 0);
        for(int vertex = 0; vertex < g->numberOfVertices; vertex++) {
            if(vertex == 0) { continue; }
            forbidEdge(g, vertex, 0);
        }
    }

    if(options->bipartiteFlag) {

        //  Since there is an n-2 cycle, we need even order.
        if(g->numberOfVertices % 2 == 1) {
            return;
        }

        //  We have already completely determined the colour classes of our
        //  bipartite graph. We now simply need to forbid edges between
        //  vertices in the same class.
        for(int evenVertex = 0; evenVertex < g->numberOfVertices - 2; evenVertex += 2) {
            for(int secondEvenVertex = evenVertex + 2; secondEvenVertex < g->numberOfVertices - 2;
             secondEvenVertex += 2) {
                forbidEdge(g, evenVertex, secondEvenVertex);
            }
            removeElement(possibilitiesForSecondEdge, evenVertex);
            forbidEdge(g, evenVertex, g->numberOfVertices - 1);
        }
        for(int oddVertex = 1; oddVertex < g->numberOfVertices - 2; oddVertex +=2) {
            for(int secondOddVertex = oddVertex + 2; secondOddVertex < g->numberOfVertices - 2;
             secondOddVertex += 2) {
                forbidEdge(g, oddVertex, secondOddVertex);
            }
            forbidEdge(g, oddVertex, g->numberOfVertices - 2);
        }
    }

    //  Dealing with girth.
    if(options->minimalGirth >= 4) {
        removeElement(possibilitiesForSecondEdge,0);
        forbidEdge(g, g->numberOfVertices - 1 , 0);
    }
    for(int k = 2; k + 4 <= options->minimalGirth; k++) {
        removeElement(possibilitiesForSecondEdge, k);
        forbidEdge(g, g->numberOfVertices - 1 , k);
        removeElement(possibilitiesForSecondEdge, g->numberOfVertices - 2 - k);  
        forbidEdge(g, g->numberOfVertices - 1 , g->numberOfVertices - 2 - k);
    }
    forEach(i, possibilitiesForSecondEdge) {
        addEdge(g, g->numberOfVertices - 1, i);

        //  Check isomorphism
        graph* gCan = malloc(sizeof(graph)*g->numberOfVertices);
        if(gCan == NULL) {
            fprintf(stderr, "Error: out of memory\n");
            exit(1);
        }
        createCanonicalForm(g->nautyGraph, gCan, g->numberOfVertices);
        bool isPresent;
        splay_insert(&splayTreeArray[g->numberOfEdges], gCan, g->numberOfVertices, &isPresent);
        if(isPresent) {
            forbidEdge(g, g->numberOfVertices - 1, i);
            free(gCan);
            removeEdge(g, g->numberOfVertices - 1, i);
            continue;
        }
        int forbiddenDueToSecondEdge[(g->numberOfVertices)*(g->numberOfVertices-1)];
        int forbiddenDueToSecondEdgeLength = 0;

        //  Third edge goes from n-2 to the cycle. It cannot join the cycle
        //  next to i. Depending on the girth restriction it cannot join
        //  the cycle girth - 4 places next to i and girth - 3 places next to 0.
        bitset possibilitiesForThirdEdge = compl(singleton(0), g->numberOfVertices - 2);
        removeElement(possibilitiesForThirdEdge, 
         (g->numberOfVertices - 2 + i - 1)%(g->numberOfVertices - 2));
        removeElement(possibilitiesForThirdEdge, (i + 1)%(g->numberOfVertices - 2));
        addForbiddenEdge(g, g->numberOfVertices - 2 ,
         (g->numberOfVertices - 2 + i - 1)%(g->numberOfVertices - 2), forbiddenDueToSecondEdge, 
         &forbiddenDueToSecondEdgeLength);
        addForbiddenEdge(g, g->numberOfVertices - 2 , (i + 1)%(g->numberOfVertices - 2),
         forbiddenDueToSecondEdge, &forbiddenDueToSecondEdgeLength);

        if(options->maximumDegree == 3) {

            //  We already do not allow 0.
            removeElement(possibilitiesForThirdEdge, i);
            for (int vertex = 0; vertex < g->numberOfVertices; ++vertex) {
                if(vertex == i) {
                    continue;
                }
                addForbiddenEdge(g, vertex, i, forbiddenDueToSecondEdge,
                 &forbiddenDueToSecondEdgeLength);
            }
        }

        if(options->bipartiteFlag) {
            for(int oddVertex = 1; oddVertex < g->numberOfVertices - 2; oddVertex +=2) {
                removeElement(possibilitiesForThirdEdge, oddVertex);
            }
        }

        for(int k = 0; k + 4 <= options->minimalGirth; k++) {
            removeElement(possibilitiesForThirdEdge, 
             (i+k)%(g->numberOfVertices - 2));
            removeElement(possibilitiesForThirdEdge, 
             (g->numberOfVertices - 2 + i - k)%(g->numberOfVertices - 2));
            removeElement(possibilitiesForThirdEdge, k + 1);
            removeElement(possibilitiesForThirdEdge,
             g->numberOfVertices - 2 - k - 1);
            addForbiddenEdge(g, g->numberOfVertices - 2 , (i+k)%(g->numberOfVertices - 2),
             forbiddenDueToSecondEdge, &forbiddenDueToSecondEdgeLength);
            addForbiddenEdge(g, g->numberOfVertices - 2 ,
             (g->numberOfVertices - 2 + i - k)%(g->numberOfVertices - 2), forbiddenDueToSecondEdge, 
             &forbiddenDueToSecondEdgeLength);
            addForbiddenEdge(g, g->numberOfVertices - 2 , k + 1, forbiddenDueToSecondEdge,
             &forbiddenDueToSecondEdgeLength);
            addForbiddenEdge(g, g->numberOfVertices - 2 , g->numberOfVertices - 2 - k - 1,
             forbiddenDueToSecondEdge, &forbiddenDueToSecondEdgeLength);
        }
        forEachAfterIndex(j, possibilitiesForThirdEdge, 0) {
            addEdge(g, g->numberOfVertices - 2, j);

            //  Check isomorphism
            graph* gCan2 = malloc(sizeof(graph)*g->numberOfVertices);
            if(gCan2 == NULL) {
                fprintf(stderr, "Error: out of memory\n");
                exit(1);
            }
            createCanonicalForm(g->nautyGraph, gCan2, g->numberOfVertices);
            bool isPresent;
            splay_insert(&splayTreeArray[g->numberOfEdges], gCan2, g->numberOfVertices, &isPresent);
            if(isPresent) {
                addForbiddenEdge(g, g->numberOfVertices - 2, j, forbiddenDueToSecondEdge,
                 &forbiddenDueToSecondEdgeLength);
                free(gCan2);
                removeEdge(g, g->numberOfVertices - 2, j);
                continue;
            }
            int forbiddenDueToThirdEdge[(g->numberOfVertices)*(g->numberOfVertices-1)];
            int forbiddenDueToThirdEdgeLength = 0;

            //  Fourth edge goes from n-2 to the cycle. It cannot join the
            //  cycle next to j or 0. Depending on the girth restriction,
            //  it cannot join the cycle girth - 4 places next to 0 and j
            //  and girth - 3 places next to i.
            bitset possibilitiesForFourthEdge = compl(singleton(i), g->numberOfVertices - 2);
            removeElement(possibilitiesForFourthEdge, g->numberOfVertices - 3);
            removeElement(possibilitiesForFourthEdge, 1);
            removeElement(possibilitiesForFourthEdge,
             (j + 1)%(g->numberOfVertices - 2));
            removeElement(possibilitiesForFourthEdge,
             (g->numberOfVertices - 2 + j - 1)%(g->numberOfVertices - 2));
            addForbiddenEdge(g, g->numberOfVertices - 1, g->numberOfVertices - 3,
             forbiddenDueToThirdEdge, &forbiddenDueToThirdEdgeLength);
            addForbiddenEdge(g, g->numberOfVertices - 1, 1, forbiddenDueToThirdEdge,
             &forbiddenDueToThirdEdgeLength);
            addForbiddenEdge(g, g->numberOfVertices - 1, (j + 1)%(g->numberOfVertices - 2),
             forbiddenDueToThirdEdge, &forbiddenDueToThirdEdgeLength);
            addForbiddenEdge(g, g->numberOfVertices - 1,
             (g->numberOfVertices - 2 + j - 1)%(g->numberOfVertices - 2), 
             forbiddenDueToThirdEdge, &forbiddenDueToThirdEdgeLength);

            if(options->maximumDegree == 3) {

                //  We already do not allow i
                removeElement(possibilitiesForFourthEdge, 0);
                removeElement(possibilitiesForThirdEdge, j);
                for (int vertex = 0; vertex < g->numberOfVertices; ++vertex)
                {
                    if(vertex == j) {
                        continue;
                    }         
                    addForbiddenEdge(g, vertex, j, forbiddenDueToThirdEdge,
                     &forbiddenDueToThirdEdgeLength);
                }
            }
            if(options->bipartiteFlag) {
                for(int evenVertex = 0; evenVertex < g->numberOfVertices - 2; evenVertex += 2) {
                    removeElement(possibilitiesForFourthEdge, evenVertex);
                }
            }

            for(int k = 0; k + 4 <= options->minimalGirth; k++) {
                removeElement(possibilitiesForFourthEdge, k);
                removeElement(possibilitiesForFourthEdge,
                 g->numberOfVertices - 2 - k);
                removeElement(possibilitiesForFourthEdge,
                 (j+k)%(g->numberOfVertices - 2));
                removeElement(possibilitiesForFourthEdge,
                 (g->numberOfVertices - 2 + j - k)%(g->numberOfVertices - 2));
                removeElement(possibilitiesForFourthEdge,
                 (i+k+1)%(g->numberOfVertices - 2));
                removeElement(possibilitiesForFourthEdge,
                 (g->numberOfVertices - 2 + i - k - 1)%(g->numberOfVertices - 2));
                addForbiddenEdge(g, g->numberOfVertices - 1, k, forbiddenDueToThirdEdge,
                 &forbiddenDueToThirdEdgeLength);
                addForbiddenEdge(g, g->numberOfVertices - 1, g->numberOfVertices - 2 - k,
                 forbiddenDueToThirdEdge, &forbiddenDueToThirdEdgeLength);
                addForbiddenEdge(g, g->numberOfVertices - 1, (j+k)%(g->numberOfVertices - 2),
                 forbiddenDueToThirdEdge, &forbiddenDueToThirdEdgeLength);
                addForbiddenEdge(g, g->numberOfVertices - 1,
                 (g->numberOfVertices - 2 + j - k)%(g->numberOfVertices - 2), 
                 forbiddenDueToThirdEdge, &forbiddenDueToThirdEdgeLength);
                addForbiddenEdge(g, g->numberOfVertices - 1, (i+k+1)%(g->numberOfVertices - 2),
                 forbiddenDueToThirdEdge, &forbiddenDueToThirdEdgeLength);
                addForbiddenEdge(g, g->numberOfVertices - 1,
                 (g->numberOfVertices - 2 + i - k - 1)%(g->numberOfVertices - 2), 
                 forbiddenDueToThirdEdge, &forbiddenDueToThirdEdgeLength);

            }
            bitset validEdges[g->numberOfVertices];
            for (int i = 0; i < (g->numberOfVertices); ++i) {
                validEdges[i] = EMPTY;
            }
            forEachAfterIndex(l, possibilitiesForFourthEdge, i) {
                if(givesValidSuccessor(g, counters, options, forbiddenDueToThirdEdge,
                 &forbiddenDueToThirdEdgeLength, validEdges, g->numberOfVertices - 1, l)) {
                    addEdge(g, g->numberOfVertices - 1, l);
                    int forbiddenDueToFourthEdge[(g->numberOfVertices)*(g->numberOfVertices-1)];
                    int forbiddenDueToFourthEdgeLength = 0;
                    addForbiddenEdge(g, g->numberOfVertices - 2, (l + 1)%(g->numberOfVertices - 2),
                     forbiddenDueToFourthEdge, &forbiddenDueToFourthEdgeLength);
                    addForbiddenEdge(g, g->numberOfVertices - 2,
                     (g->numberOfVertices - 2 + l - 1)%(g->numberOfVertices - 2), 
                     forbiddenDueToFourthEdge, &forbiddenDueToFourthEdgeLength);

                    if(options->maximumDegree == 3) {
                        for(int vertex = 0; vertex < g->numberOfVertices; vertex++) {
                            if(vertex == l) {
                                continue;
                            }
                            addForbiddenEdge(g, vertex, l, forbiddenDueToFourthEdge,
                             &forbiddenDueToFourthEdgeLength);
                        }
                    }

                    for(int k = 0; k + 4 <= options->minimalGirth; k++) {
                        addForbiddenEdge(g, g->numberOfVertices - 2,
                         (l+k)%(g->numberOfVertices - 2), forbiddenDueToFourthEdge, 
                         &forbiddenDueToFourthEdgeLength);
                        addForbiddenEdge(g, g->numberOfVertices - 2,
                         (g->numberOfVertices - 2 + l - k)%(g->numberOfVertices - 2), 
                         forbiddenDueToFourthEdge, &forbiddenDueToFourthEdgeLength);
                    }

                    int newForbiddenEdges[g->numberOfVertices*(g->numberOfVertices - 1)];
                    int newForbiddenEdgesLength = 0;

                    // writeToG6(g->nautyGraph, g->numberOfVertices);
                    // printGraph(g);

                    extendWithEdge(g, counters, options,
                     g->numberOfVertices - 1, l, newForbiddenEdges,
                     &newForbiddenEdgesLength);
                    permitEdges(g, newForbiddenEdges, newForbiddenEdgesLength);
                    permitEdges(g, forbiddenDueToFourthEdge, forbiddenDueToFourthEdgeLength);
                    removeEdge(g, g->numberOfVertices - 1, l);
                }
            }
            removeEdge(g, g->numberOfVertices - 2, j);

            //  Permit edges forbidden due to (n-2, j).
            permitEdges(g, forbiddenDueToThirdEdge, forbiddenDueToThirdEdgeLength);
        }
        removeEdge(g, g->numberOfVertices - 1, i);

        //  Permit the edges which were forbidden due to (i,n-1).
        permitEdges(g, forbiddenDueToSecondEdge, forbiddenDueToSecondEdgeLength);
    }
}

//  Generate K2hypohamiltonian graphs of order numberOfVertices.
void generateK2HypohamiltonianGraphs(int numberOfVertices, 
    struct counters *counters, struct options *options) {
    struct graph g = {.numberOfVertices = numberOfVertices};
    g.adjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    g.forbiddenEdges = malloc(sizeof(bitset)*numberOfVertices);
    g.verticesOfDeg = malloc(sizeof(bitset)*numberOfVertices);
    for(int i = 0; i < g.numberOfVertices; i++) {
        g.verticesOfDeg[i] = EMPTY;
    }
    emptyGraph(&g);

    // Create (n-2)-cycle.
    for (int i = 1; i < numberOfVertices - 2; i++) {
        addEdge(&g, i - 1, i);
    }
    addEdge(&g, numberOfVertices - 3, 0);

    // Add isolated K2.
    addEdge(&g, numberOfVertices - 2, numberOfVertices - 1);

    addK2Edges(&g, counters, options);
    free(g.adjacencyList);
    free(g.forbiddenEdges);
    free(g.verticesOfDeg);
}

int main(int argc, char ** argv) {

    /*Initialize the optional command line flags.*/

    int opt;
    bool haveModResPair = false;
    bool printCountsPerEdgeNumber = false;
    struct options options = 
        {.minimalGirth = -1,
         .minimumDegree = -1,
         .maximumDegree = -1,
         .planarFlag = false,
         .bipartiteFlag = false,
         .splitLevel = SPLITLEVEL,
         .splitCounter = 0,
         .remainder = 0,
         .modulo = 1};
    char *graphClass = "";
    while(1) {
        int option_index = 0;
        static struct option long_options[] =
        {
            {"bipartite", no_argument, NULL, 'b'},
            {"minimum-degree", required_argument, NULL, 'd'},
            {"maximum-degree", required_argument, NULL, 'D'},
            {"edges-count", no_argument, NULL, 'e'},
            {"girth", required_argument, NULL, 'g'},
            {"help", no_argument, NULL, 'h'},
            {"planar", no_argument, NULL, 'p'},
            {"splitlevel", required_argument, NULL, 'X'}
        };
        opt = getopt_long(argc, argv, "bcd:D:eg:hpX:", long_options, &option_index);
        if(opt == -1) break;
        char *ptr;
        switch(opt) {
            case 'b':
                options.bipartiteFlag = true;
                graphClass = "bipartite ";
                break;
            case 'd':
                options.minimumDegree = strtol(optarg, &ptr, 10);
                if(options.minimumDegree < 3) {
                    fprintf(stderr, 
                     "Warning: the minimal degree of a K2-hypohamiltonian graph is at least 3.\n");
                } 
                break;
            case 'D':
                options.maximumDegree = strtol(optarg, &ptr, 10);
                if(options.maximumDegree < 3) {
                    fprintf(stderr,
                     "Error: K2-hypohamiltonian graphs have degree at least 3.\n");
                    return 1;
                }
                break;
            case 'e':
                printCountsPerEdgeNumber = true;
                break;
            case 'h':
                fprintf(stderr, "%s", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'g': 
                options.minimalGirth = strtol(optarg, &ptr, 10);
                if(options.minimalGirth < 3) {
                    fprintf(stderr, "Error: the minimal girth needs to be at least 3.\n");
                    fprintf(stderr, "%s", USAGE);
                    return 1;
                }
                break;
            case 'p':
                options.planarFlag = true;
                graphClass = "planar ";
                break;
            case 'X':
                options.splitLevel = strtol(optarg, &ptr, 10);
                break;
        }
    }
    //  Check if there is a non-option argument.
    if (optind >= argc) {
        fprintf(stderr, "Error: add number of vertices.");
        fprintf(stderr, "%s", USAGE);
        return 1;
    }
    //  First non-option argument should be number of vertices.
    char* endptr;
    int numberOfVertices = strtol(argv[optind++], &endptr, 10);

    if(numberOfVertices <= 2 || numberOfVertices > MAXN) {
        fprintf(stderr, "Error: n needs to be a number between 3 and %d.",
         MAXN);
        fprintf(stderr, "%s",USAGE);
        return 1;
    }

    //  Check for other non-option arguments.
    while (optind < argc) {
        bool pairIsInvalid = false;
        options.remainder = strtol(argv[optind], &endptr, 10);
        if( !endptr || *endptr != '/' || *(endptr+1) == '\0') {
            pairIsInvalid = true;
        }
        options.modulo = strtol(endptr+1, &endptr, 10);
        if( !endptr || *endptr != '\0') {
            pairIsInvalid = true;
        }
        if(options.modulo <= options.remainder) {
            pairIsInvalid = true;
        }
        if(haveModResPair) {
            fprintf(stderr,
             "Error: You can only add one res/mod pair as an argument.\n");
            fprintf(stderr, "%s\n", USAGE);
            fprintf(stderr,
             "Use ./genK2Hypohamiltonian --help for more detailed instructions.\n");
            return 1;
        }
        haveModResPair = true;

        if(pairIsInvalid) {
            fprintf(stderr,
                 "Error: Invalid res/mod pair: '%s'.\n", argv[optind]);
            fprintf(stderr, "%s\n", USAGE);
            fprintf(stderr,
             "Use ./genK2Hypohamiltonian --help for more detailed instructions.\n");
            return 1;
        }
        fprintf(stderr, "Class=%d/%d. Splitlevel = %d.\n",
         options.remainder, options.modulo, options.splitLevel);
        optind++;
    }

    /*Do some correctness checks.*/

    //  Cubic K2-hypohamiltonian graphs have girth at least 5.
    if(options.maximumDegree == 3) {
        if(options.minimalGirth < 5) {
            options.minimalGirth = 5;
        }
    }

    //  Bipartite K2-hypohamiltonian graphs can only have even girth.
    if(options.bipartiteFlag && options.minimalGirth != -1) {
        options.minimalGirth += options.minimalGirth%2;
    }
    nauty_check(WORDSIZE, SETWORDSNEEDED(numberOfVertices),
     numberOfVertices, NAUTYVERSIONID);

    /*Initialize all counters to 0.*/
    struct counters counters = {};
    counters.nOfNonIsoGraphsWithEdges = calloc(sizeof(long long unsigned int),
     numberOfVertices*(numberOfVertices - 1)/2+1);

    /*Main routine*/
    clock_t start = clock();
    generateK2HypohamiltonianGraphs(numberOfVertices, &counters, &options);
    clock_t end = clock();

    /*Output results*/
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    if(options.minimalGirth != -1) {
        fprintf(stderr,
         "\rGenerated %lld %sgraphs of order %d with girth at least %d in %f seconds.\n",
         counters.nOfGraphsFound, graphClass, numberOfVertices, options.minimalGirth, time_spent);
    }
    else {
        fprintf(stderr,
         "\rGenerated %lld %sgraphs of order %d in %f seconds.\n", 
         counters.nOfGraphsFound, graphClass, numberOfVertices, time_spent);
    }
    if(options.minimumDegree != -1) {
        fprintf(stderr, "Their minimum degree is %d.\n", options.minimumDegree);
    }
    if(options.maximumDegree != -1) {
        fprintf(stderr, "Their maximum degree is %d.\n", options.maximumDegree);
    }
    if(options.planarFlag) {
        fprintf(stderr, "Times graph was non-planar: %llu (%.2f%%)\n",
         counters.nOfTimesWasNonPlanar, 
         (double) 100*counters.nOfTimesWasNonPlanar/counters.nOfTimesCheckedPlanarity);
    }
    fprintf(stderr, "---\n");
    fprintf(stderr, "Number of non-isomorphic graphs checked: %llu\n",
     counters.nOfTimesCheckedIsomorphism - counters.nOfTimesWasIsomorphic);
    fprintf(stderr, "Times contained type A: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n",
     counters.nOfTimesContainedTypeA, 
     (double) 100*counters.nOfTimesContainedTypeA/counters.nOfTimesCheckedTypeA,
     counters.nOfTimesTypeAObstructionChosen, 
     (double) 100*counters.nOfTimesTypeAObstructionChosen/counters.nOfTimesObstructionFound);
    // fprintf(stderr, "Times contained type B: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n",
    //  counters.nOfTimesContainedTypeB, 
    //  (double) 100*counters.nOfTimesContainedTypeB/counters.nOfTimesCheckedTypeB,
    //  counters.nOfTimesTypeBObstructionChosen, 
    //  (double) 100*counters.nOfTimesTypeBObstructionChosen/counters.nOfTimesObstructionFound);
    fprintf(stderr, "Times contained type C: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n",
     counters.nOfTimesContainedTypeC, 
     (double) 100*counters.nOfTimesContainedTypeC/counters.nOfTimesCheckedTypeC,
     counters.nOfTimesTypeCObstructionChosen, 
     (double) 100*counters.nOfTimesTypeCObstructionChosen/counters.nOfTimesObstructionFound);
    fprintf(stderr, "Times contained bad degree vertex: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n",
     counters.nOfTimesContainedDegreeObstruction, 
     (double) 100*counters.nOfTimesContainedDegreeObstruction/counters.nOfTimesCheckedDegreeObstruction,
     counters.nOfTimesDegreeObstructionChosen, 
     (double) 100*counters.nOfTimesDegreeObstructionChosen/counters.nOfTimesObstructionFound);
    if(options.minimalGirth < 3) {
        fprintf(stderr,
         "Times contained general triangle obstruction: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n", 
         counters.nOfTimesContainedStarObstruction, 
         (double) 100*counters.nOfTimesContainedStarObstruction/counters.nOfTimesCheckedStarObstruction,
         counters.nOfTimesGeneralTriangleObstructionChosen, 
         (double) 100*counters.nOfTimesGeneralTriangleObstructionChosen/counters.nOfTimesObstructionFound);   
        fprintf(stderr,
         "Times contained deg 3 triangle vertex: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n", 
         counters.nOfTimesContainedTriangleObstruction, 
         (double) 100*counters.nOfTimesContainedTriangleObstruction/counters.nOfTimesCheckedTriangleObstruction,
         counters.nOfTimesTriangleObstructionChosen, 
         (double) 100*counters.nOfTimesTriangleObstructionChosen/counters.nOfTimesObstructionFound);
    }
    if(options.minimalGirth < 4) {
        fprintf(stderr,
         "Times contained deg 3 4-cycle vertex: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n", 
         counters.nOfTimesContainedSquareObstruction, 
         (double) 100*counters.nOfTimesContainedSquareObstruction/counters.nOfTimesCheckedSquareObstruction,
         counters.nOfTimes4CycleObstructionChosen, 
         (double) 100*counters.nOfTimes4CycleObstructionChosen/counters.nOfTimesObstructionFound);
        fprintf(stderr,
         "Times contained general 4-cycle obstruction: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n", 
         counters.nOfTimesContainedArrowObstruction, 
         (double) 100*counters.nOfTimesContainedArrowObstruction/counters.nOfTimesCheckedArrowObstruction,
         counters.nOfTimesGeneral4CycleObstructionChosen, 
         (double) 100*counters.nOfTimesGeneral4CycleObstructionChosen/counters.nOfTimesObstructionFound); 
    }   
    fprintf(stderr, "Times no obstruction applied: %llu (%.2f%% of non-isomorphic graphs)\n",
     counters.nOfTimesNoObstructionApplied, 
     (double) 100*counters.nOfTimesNoObstructionApplied / 
     (counters.nOfTimesCheckedIsomorphism - counters.nOfTimesWasIsomorphic));
    fprintf(stderr, "---\n");
    fprintf(stderr, "Times adding successor edge made graph hamiltonian: %llu (%.2f%%)\n",
     counters.nOfTimesWasHamiltonian, 
     (double) 100*counters.nOfTimesWasHamiltonian/counters.nOfTimesCheckedHamiltonicity);
    if(options.minimalGirth >= 4) {
        fprintf(stderr, "Times adding successor edge gave forbidden girth: %llu (%.2f%%)\n",
         counters.nOfTimesHadForbiddenGirth, 
         (double) 100*counters.nOfTimesHadForbiddenGirth/counters.nOfTimesCheckedGirth);
    }
    if(options.maximumDegree != -1) {
        fprintf(stderr, "Times adding successor edge gave vertex forbidden degree: %llu (%.2f%%)\n",
         counters.nOfTimesMaximumDegreeExceeded, 
         (double) 100*counters.nOfTimesMaximumDegreeExceeded/counters.nOfTimesCheckedMaximumDegree);
    }
    fprintf(stderr, "Size of graphs in splay trees: %.3f gigabytes\n", 
     (double) (sizeof(graph)*numberOfVertices *
     (counters.nOfTimesCheckedIsomorphism - counters.nOfTimesWasIsomorphic))/1000000000);

    if(printCountsPerEdgeNumber) {
        for(int i = numberOfVertices+3; i <= numberOfVertices*(numberOfVertices - 1)/2; i++) {
            fprintf(stderr, "%8llu non-isomorphic graphs with %d edges generated.\n",
             counters.nOfNonIsoGraphsWithEdges[i], i);
            if(counters.nOfNonIsoGraphsWithEdges[i]==0) {
                fprintf(stderr,"No graphs with more edges generated.\n");
                break;
            }
        }
    }
    return 0;
}
