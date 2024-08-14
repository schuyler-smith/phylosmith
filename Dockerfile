FROM r-base:4.4.1

WORKDIR /opt
ENV PATH="${PATH}:/opt/"

RUN apt-get update \
  && apt-get install --yes --no-install-recommends \
    ca-certificates \
    build-essential \
    zip unzip\
    wget \
    gfortran \
    g++ \
    libbz2-dev \
    libpq-dev \
    libgdal-dev \
    libudunits2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libssl-dev \
    libxml2-dev \
  && apt-get upgrade --yes \
  && apt-get clean

  RUN apt-get update \
  && apt-get install --yes --no-install-recommends --fix-broken --allow-downgrades \
    libcurl4-openssl-dev \
    libcurl4t64=8.8.0-4 \
  && apt-get upgrade --yes \
  && apt-get clean

RUN R -e "install.packages('remotes', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('knitr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('biomformat')"
RUN R -e "BiocManager::install('Biostrings')"

COPY DESCRIPTION DESCRIPTION
RUN Rscript -e 'desc <- read.dcf("DESCRIPTION"); deps <- c(desc[1, "Imports"], desc[1, "LinkingTo"]); install.packages(trimws(unlist(strsplit(deps, ","))), repos="http://cran.rstudio.com/")'


WORKDIR /app
ENV PATH="${PATH}:/app/"

RUN R -e "BiocManager::install('phyloseq')"
COPY . .
RUN R CMD INSTALL .

ENTRYPOINT ["/bin/bash"]

