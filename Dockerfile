## Parent image
FROM rocker/r-ubuntu:20.04


## Make and set directories
RUN mkdir /opt/genepatt
WORKDIR /opt/genepatt
RUN mkdir src
RUN mkdir testdata

## make directory for packages, and copy them in
RUN mkdir pkgs
COPY packages/* pkgs/



## install packages
RUN Rscript -e "install.packages('optparse', version = '1.7.3', repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('BiocManager', version = '1.30.18')"
RUN Rscript -e "library('BiocManager')"
### edgeR version: 3.38.4
RUN Rscript -e "BiocManager::install('edgeR')"


## copying data
COPY data/* /opt/genepatt/testdata/

## Copying source code
COPY src/* /opt/genepatt/src/
