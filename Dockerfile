## Parent image
FROM r-base:4.2.2


## Make and set directories
RUN mkdir /opt/genepatt
WORKDIR /opt/genepatt
RUN mkdir src
RUN mkdir testdata



## install packages
RUN Rscript -e "install.packages('optparse', version = '1.7.3', repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('BiocManager', version = '1.30.18')"
RUN Rscript -e "library('BiocManager')"
### edgeR version: 3.38.4
RUN Rscript -e "BiocManager::install('edgeR')"


## install python
RUN apt update
RUN apt install -y python3

## copying data
COPY data/* /opt/genepatt/testdata/

## Copying source code
COPY src/* /opt/genepatt/src/
