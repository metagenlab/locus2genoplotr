

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
singularity build locus2genoplotr.simg docker://metagenlab/locus2genoplotr:1.0
```

### example

```bash
singularity exec locus2genoplotr.py -l KL28_00008 -r KL28.gbk -q capsule_region_150bp_assembly_concat.gbk capsule_region_250bp_assembly_concat.gbk -rs 45000 -ls 30000
```