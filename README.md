

# Description

Plot gene or genome maps and comparisons genome plots. Generate *SVG* and *PNG* images with genoPlotR (http://genoplotr.r-forge.r-project.org/) or GenomeDiagram (https://academic.oup.com/bioinformatics/article/22/5/616/205776).

# Installation

## Method 1: Installation with conda

```bash
# create the environment
conda env create -f env.yaml
# activate the environment
conda activate locus2genoplotr
```

## Method 2: Singularity container   

Dependency: [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)

### build the image 

```bash
singularity build locus2genoplotr.simg docker://metagenlab/locus2genoplotr:1.1
```

## examples


### simple plot 

![Simple plot](examples/plasmid_linear.svg)

### comparative plot

```bash
locus2genoplotr -l KL28_00008 -r data/KL28.gbk -q data/capsule_region_150bp_assembly_concat.gbk data/capsule_region_250bp_assembly_concat.gbk -rs 45000 -ls 30000 -o simple_comp -v
```
![Comparison plot](examples/simple_comp.svg)

### comparative plot tblastx

```bash
locus2genoplotr.py -l KL24_00002 -r data/K24.gbk -q data/GCF_001596925.1_ASM159692v1_genomic.gbff data/GCF_000943095.1_ST15_genomic.gbff data/K24.gbk -rs 25000 -ls 3000 -x -o alignment_tblastx
```

![TblastX plot](examples/alignment_tblastx.svg)

### comparative plot with GC 

```bash
locus2genoplotr -l KL28_00008 -r data/KL28.gbk -q data/capsule_region_150bp_assembly_concat.gbk data/capsule_region_250bp_assembly_concat.gbk -rs 45000 -ls 30000 -s data/capsule_region_250bp_assembly_concat.depth_150bp_reads -o capsule_with_gc -v
```

![GC plot](examples/capsule_with_gc.svg)

### comparative plot with sequencing depth 

```bash
locus2genoplotr -l KL28_00008 -r data/KL28.gbk -q data/capsule_region_150bp_assembly_concat.gbk data/capsule_region_250bp_assembly_concat.gbk -rs 45000 -ls 30000 -o capsule_with_gc -v -g
```

![Depth plot](examples/capsule_with_depth.svg)