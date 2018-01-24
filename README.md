

# installation

## python

```
# set conda priorities
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

conda install -c conda-forge readline=6.2
#conda install readline
conda install rpy2
```
## R

```
conda install r-base
conda install r-ggplot2
```

```
install.packages("genoPlotR")
```
