# kmotif

## Info:

Efficient C code related to listing k-motifs in large sparse graphs. 
- k-motif: connected induced subgraph on k nodes.
- k-path: k-motif such that k-2 nodes have degree 2 and 2 nodes have degree 1.
- k-star: k-motif such that one of the nodes is adjacent to all other k-1 nodes (it does not have to be a tree). Maybe I'll find another name for this...

Feel free to use these lines as you wish.

## To compile:

"gcc kmotif.c -O3 -o kmotif"  
"gcc kpath.c -O3 -o kpath"  
"gcc rmhub.c -O3 -o rmhub"  
"gcc kstarLBUB -O3 -o kstarLBUB"

## To execute:

"./kmotif k edgelist.txt".
- "k" is the size (number of nodes) of the k-motifs to consider
- "edgelist.txt" should contain the graph: one edge on each line separated by a space (no self-loops and no multiple edges).
- Will print the number of k-motifs.

"./kpath k edgelist.txt".
- "k" is the size (number of nodes) of the k-paths to consider
- "edgelist.txt" should contain the graph: one edge on each line separated by a space (no self-loops and no multiple edges).
- Will print the number of k-paths.

Hubs can be removed (if the program is too slow) using the following command.  
"./rmhub dmax input.txt output.txt".
- "dmax" should be an integer (maximum degree)
- "input.txt" should contain the graph: one edge on each line separated by a space.
- "output.txt" is the graph with all nodes of degree more than dmax removed.

"./kstarLBUB k edgelist.txt".
- "k" is the size (number of nodes) of the k-motifs to consider
- "edgelist.txt" should contain the graph: one edge on each line separated by a space (no self-loops and no multiple edges).Will - Will print an lower-bound and upper-bound on the number of k-stars.



## Performance (on a commodity machine):

On https://snap.stanford.edu/data/com-Amazon.html:
- 8,496,678,114 5-motifs in 9 seconds; 250,355,118 5-paths in 7 seconds
- 623,221,614,256 6-motifs in 4 minutes; 1,651,951,945 6-paths in 1 minute
- 11,804,852,312 7-paths in 4 minutes

On https://snap.stanford.edu/data/com-Youtube.html:
- 5,834,725,962,171 4-motifs in 3 minutes; 91,488,735,459 4-paths in 12 minutes

On https://snap.stanford.edu/data/com-Orkut.html:
- 44,370,298,209 3-motifs in 3 minutes; 43,742,714,028 3-paths in 6 minutes (the number of triangles in Orkut is indeed the difference: 627,584,181).

On https://snap.stanford.edu/data/com-Friendster.html:
- 712,307,515,908 3-motifs in 2 hours; 708,133,791,784 3-paths in 5 hours (the number of triangles in Friendster is indeed the difference: 4,173,724,124)
- 618,294,018,778,462 4-motifs in 610 hours;

Note that for motifs the reported time is the time to count the motifs and not the time to list them which is a bit longueur (c.f. line 205-225 of kmotif.c). For paths, the reported time is the time to list all paths...  
All codes can be made parallel.


## Initial contributors:

Maximilien Danisch  
March 2018  
http://bit.ly/danisch  
maximilien.danisch@gmail.com
