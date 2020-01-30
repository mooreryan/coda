# coda

Do a couple of essential compositional data analysis things on giant data sets.

## Build

### Prereqs

* CMake 3.15
* C++17

### Install from source

Download the source code, then run 

```
mkdir build
cd build
cmake ..
make
make install
```

## Running

This program assumes that OTUs (or whatever) are on rows, and samples are on columns.  Also, it assumes that there are fewer samples than OTUs/clusters/seqs.

*Note: `seed` is the random seed used for Randomized SVD.*

```
coda <seed> <counts> <clr_out.tsv> <aitchison_dist_out.tsv> <sample_projection_out.tsv>
```

If you have OpenMP, then you can set number of threads like this:

```
OMP_NUM_THREADS=n coda ...
```
 
Where `coda ...` is the command you would normally run.

## Output 

* CLR transformed OTU table
* Aitchison distance matrix for the samples
* Projection of samples into PC space based on the CLR transformation (like doing a compositional biplot)
