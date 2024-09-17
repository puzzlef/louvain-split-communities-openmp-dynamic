Design of OpenMP-based Parallel [Louvain algorithm][Louvain] for [community detection], \
with no internally-disconnected communities.

<br>

Community detection entails the identification of clusters of vertices that exhibit stronger connections within themselves compared to the wider network. The Louvain method, a commonly utilized heuristic for this task, employs a two-step process comprising a local-moving phase and an aggregation phase. This process iteratively optimizes the modularity metric, a measure of community quality. Although widely used, the Louvain method has been noted for generating internally-disconnected and poorly connected communities. In response, Traag et al. introduce the Leiden algorithm, which incorporates an extra refinement phase between the local-moving and aggregation phases. This refinement step enables vertices to explore and potentially create sub-communities within the identified communities from the local-moving phase. In this report, we propose another BFS-based approach, **GSP-Louvain**, for addressing the same issue. This is different from a number of earlier works, which tackled internally-disconnected communities as a post-processing step.

Below we plot the time taken by the [original Leiden], [igraph] Leiden, [NetworKit] Leiden, GSP-Louvain on 13 different graphs. GSP-Louvain surpasses the original Leiden, igraph Leiden, and NetworKit Leiden by `341×`, `83×`, and `6.1×` respectively, achieving a processing rate of `328M` edges/s on a `3.8𝐵` edge graph.

[![](https://i.imgur.com/Ed9Slnw.png)][sheets-o1]

Below we plot the speedup of GSP-Louvain wrt original Leiden, igraph Leiden, and NetworKit Leiden.

[![](https://i.imgur.com/0M99mgH.png)][sheets-o1]

Next, we compare the modularity of communities identified by the original Leiden algorithm, igraph Leiden, NetworKit Leiden, and GSP-Louvain. On average, GSP-Louvain achieves `0.3%` lower modularity than the original Leiden and igraph Leiden, respectively, and `25%` higher modularity than NetworKit Leiden, particularly evident on road networks and protein k-mer graphs.

[![](https://i.imgur.com/GT6vFxZ.png)][sheets-o1]

Finally, we plot the fraction of disconnected communities identified by each implementation. Absence of bars indicates the absence of disconnected communities. The original Leiden, igraph Leiden, and GSP-Louvain do not identify any communities that are internally-disconnected. However, on average, NetworKit Leiden exhibit fraction of disconnected communities amounting to `1.5×10^−2`, particularly on web graphs and social networks.

[![](https://i.imgur.com/dIDYXhP.png)][sheets-o1]

Refer to our technical reports for more details: \
[An Approach for Addressing Internally-Disconnected Communities in Louvain Algorithm][report].

<br>

> [!NOTE]
> You can just copy `main.sh` to your system and run it. \
> For the code, refer to `main.cxx`.

[Leiden]: https://www.nature.com/articles/s41598-019-41695-z
[Louvain]: https://arxiv.org/abs/0803.0476
[original Leiden]: https://github.com/vtraag/libleidenalg
[igraph]: https://github.com/igraph/igraph
[NetworKit]: https://github.com/networkit/networkit
[community detection]: https://en.wikipedia.org/wiki/Community_search
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[sheets-o1]: https://docs.google.com/spreadsheets/d/1N8eoVV5AUFYUKuvZBbHvL1BPc86xgmgAPA_t4pIC1gk/edit?usp=sharing
[report]: https://arxiv.org/abs/2402.11454

<br>
<br>


### Code structure

The code structure of GSP-Louvain is as follows:

```bash
- inc/_algorithm.hxx: Algorithm utility functions
- inc/_bitset.hxx: Bitset manipulation functions
- inc/_cmath.hxx: Math functions
- inc/_ctypes.hxx: Data type utility functions
- inc/_cuda.hxx: CUDA utility functions
- inc/_debug.hxx: Debugging macros (LOG, ASSERT, ...)
- inc/_iostream.hxx: Input/output stream functions
- inc/_iterator.hxx: Iterator utility functions
- inc/_main.hxx: Main program header
- inc/_mpi.hxx: MPI (Message Passing Interface) utility functions
- inc/_openmp.hxx: OpenMP utility functions
- inc/_queue.hxx: Queue utility functions
- inc/_random.hxx: Random number generation functions
- inc/_string.hxx: String utility functions
- inc/_utility.hxx: Runtime measurement functions
- inc/_vector.hxx: Vector utility functions
- inc/batch.hxx: Batch update generation functions
- inc/bfs.hxx: Breadth-first search algorithms
- inc/csr.hxx: Compressed Sparse Row (CSR) data structure functions
- inc/dfs.hxx: Depth-first search algorithms
- inc/duplicate.hxx: Graph duplicating functions
- inc/Graph.hxx: Graph data structure functions
- inc/louvian.hxx: Louvian algorithm functions
- inc/louvainSplit.hxx: Louvain with no disconnected communities
- inc/main.hxx: Main header
- inc/mtx.hxx: Graph file reading functions
- inc/properties.hxx: Graph Property functions
- inc/selfLoop.hxx: Graph Self-looping functions
- inc/symmetrize.hxx: Graph Symmetricization functions
- inc/transpose.hxx: Graph transpose functions
- inc/update.hxx: Update functions
- main.cxx: Experimentation code
- process.js: Node.js script for processing output logs
```

Note that each branch in this repository contains code for a specific experiment. The `main` branch contains code for the final experiment. If the intention of a branch in unclear, or if you have comments on our technical report, feel free to open an issue.

<br>
<br>


## References

- [Fast unfolding of communities in large networks; Vincent D. Blondel et al. (2008)](https://arxiv.org/abs/0803.0476)
- [Community Detection on the GPU; Md. Naim et al. (2017)](https://arxiv.org/abs/1305.2006)
- [Scalable Static and Dynamic Community Detection Using Grappolo; Mahantesh Halappanavar et al. (2017)](https://ieeexplore.ieee.org/document/8091047)
- [From Louvain to Leiden: guaranteeing well-connected communities; V.A. Traag et al. (2019)](https://www.nature.com/articles/s41598-019-41695-z)
- [CS224W: Machine Learning with Graphs | Louvain Algorithm; Jure Leskovec (2021)](https://www.youtube.com/watch?v=0zuiLBOIcsw)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [Fetch-and-add using OpenMP atomic operations](https://stackoverflow.com/a/7918281/1413259)

<br>
<br>


[![](https://i.imgur.com/atJbkL1.png)](https://www.youtube.com/watch?v=yqO7wVBTuLw&pp)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/519156419.svg)](https://zenodo.org/doi/10.5281/zenodo.6945748)


[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
