FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
ENV R_VERSION=4.4.2

WORKDIR /opt
ENV PATH="${PATH}:/opt/"

RUN apt-get update -qq \
  && apt-get install --yes --no-install-recommends \
    wget curl \
    make cmake \
    autoconf automake \
    zip unzip gzip \
    build-essential \
    gfortran \
    libreadline-dev \
    xorg-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    liblzma-dev \
    libpcre2-dev \ 
  && apt-get upgrade --yes \
  && apt-get clean

RUN apt-get install --no-install-recommends --yes \
    locales \
  && apt-get upgrade --yes; apt-get clean \
  && echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8
  ENV LC_ALL=en_US.UTF-8
  ENV LANG=en_US.UTF-8

RUN wget -c --no-check-certificate https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz \
  && tar -xf R-${R_VERSION}.tar.gz \
  && cd R-${R_VERSION} \
  && ./configure \
  && make -j$(nproc) \
  && make install \
  && cd .. \
  && rm -rf R-${R_VERSION} R-${R_VERSION}.tar.gz

# R Package Install Functions
RUN echo "#!/usr/bin/env Rscript" > /usr/local/bin/import_r_cl_packages \
  && echo "import_r_cl_packages <- function(args) {" >> /usr/local/bin/import_r_cl_packages \
  && echo "args <- gsub(\"'\", \"\", args)" >> /usr/local/bin/import_r_cl_packages \
  && echo "args <- gsub(\" \", \"\", args)" >> /usr/local/bin/import_r_cl_packages \
  && echo "packages <- unlist(strsplit(args, \",\"))" >> /usr/local/bin/import_r_cl_packages \
  && echo "return(packages)}" >> /usr/local/bin/import_r_cl_packages

RUN echo "#!/usr/bin/env Rscript" > /usr/local/bin/check_r_packages \
  && echo "source('/usr/local/bin/import_r_cl_packages')" >> /usr/local/bin/check_r_packages \
  && echo "check_r_packages <- function(x) {missing <- setdiff(x, rownames(installed.packages()))" >> /usr/local/bin/check_r_packages \
  && echo "  if (length(missing) > 0) stop('Missing: ', paste(missing, collapse=', '))" >> /usr/local/bin/check_r_packages \
  && echo "  else cat(paste(x, 'installed successfully\n'))}" >> /usr/local/bin/check_r_packages \
  && echo "packages = import_r_cl_packages(commandArgs(trailingOnly=TRUE)[1])" >> /usr/local/bin/check_r_packages \
  && echo "check_r_packages(packages)" >> /usr/local/bin/check_r_packages \
  && chmod +x /usr/local/bin/check_r_packages

RUN echo "#!/usr/bin/env Rscript" > /usr/local/bin/install_r_packages \
  && echo "source('/usr/local/bin/import_r_cl_packages')" >> /usr/local/bin/install_r_packages \
  && echo "packages = import_r_cl_packages(commandArgs(trailingOnly=TRUE)[1])" >> /usr/local/bin/install_r_packages \
  && echo "install.packages(packages, repos='http://cran.rstudio.com/')" >> /usr/local/bin/install_r_packages \
  && chmod +x /usr/local/bin/install_r_packages

RUN apt-get install --no-install-recommends --yes \
    git ca-certificates \
  && apt-get upgrade --yes; apt-get clean
RUN install_r_packages 'remotes' \
  && check_r_packages 'remotes'

ENV R_PACKAGES="'ggh4x', 'ggpubr', 'ggrepel', 'ggforce', 'Rtsne', 'RColorBrewer', 'viridis', 'viridisLite', 'igraph', 'dendextend', 'ggraph', 'data.table', 'knitr'"
RUN install_r_packages "${R_PACKAGES}" \
  && check_r_packages "${R_PACKAGES}"

RUN apt-get install --no-install-recommends --yes \
    libudunits2-dev libgdal-dev \
  && apt-get upgrade --yes; apt-get clean
RUN install_r_packages 'sf' \
  && check_r_packages 'sf'


RUN echo | openssl s_client -showcerts -servername bioconductor.org -connect bioconductor.org:443 > bioconductor.pem \
  && cp bioconductor.pem /usr/local/share/ca-certificates/bioconductor.crt \
  && update-ca-certificates
RUN install_r_packages 'BiocManager' \
  && check_r_packages 'BiocManager'

ENV R_PACKAGES="'BiocGenerics', 'Biobase', 'S4Vectors', 'ShortRead', 'IRanges', 'zlibbioc', \
  'XVector', 'UCSC.utils', 'GenomeInfoDbData', 'GenomeInfoDb', 'Biostrings', 'multtest'"
RUN R -e "BiocManager::install(c(${R_PACKAGES}))" \
  && check_r_packages "${R_PACKAGES}"

ENV R_PACKAGES="'phyloseq'"
RUN R -e "BiocManager::install(c(${R_PACKAGES}))" \
  && check_r_packages "${R_PACKAGES}"

ENV R_PACKAGES="'RCurl', 'RcppProgress', 'RcppParallel'"
RUN install_r_packages "${R_PACKAGES}" \
  && check_r_packages "${R_PACKAGES}"

# COPY DESCRIPTION DESCRIPTION
# RUN Rscript -e 'desc <- read.dcf("DESCRIPTION"); deps <- c(desc[1, "Imports"], desc[1, "LinkingTo"]); install.packages(trimws(unlist(strsplit(deps, ","))), repos="http://cran.rstudio.com/")'


RUN apt-get install --no-install-recommends --yes \
    libharfbuzz-dev libfribidi-dev \
  && apt-get upgrade --yes; apt-get clean
ENV R_PACKAGES="'devtools'"
RUN install_r_packages "${R_PACKAGES}" \
  && check_r_packages "${R_PACKAGES}"


WORKDIR /app
ENV PATH="${PATH}:/app/"

# COPY . .
# RUN R CMD INSTALL .

ENTRYPOINT ["/bin/bash"]

