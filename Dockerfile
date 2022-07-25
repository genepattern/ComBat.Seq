## Parent image
FROM rocker/r-ubuntu:20.04


## Make and set directories
RUN mkdir /opt/genepatt
WORKDIR /opt/genepatt
RUN mkdir src
RUN mkdir testdata


## Copying source code
COPY src/* /opt/genepatt/src/



## install packages


## optional packages
RUN R -e "install.packages('optparse', repos='http://cran.us.r-project.org')"


## copying data
COPY data/* /opt/genepatt/testdata/
