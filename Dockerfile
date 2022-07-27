## Parent image
FROM rocker/r-ubuntu:20.04


## Make and set directories
RUN mkdir /opt/genepatt
WORKDIR /opt/genepatt
RUN mkdir src
RUN mkdir testdata

## install packages
RUN Rscript -e "install.packages('optparse', repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "library('BiocManager')"
RUN Rscript -e "BiocManager::install(c('edgeR', 'sva'))"

## copying data
COPY data/* /opt/genepatt/testdata/

## Copying source code
COPY src/* /opt/genepatt/src/
