FROM rocker/binder:3.6.0
LABEL maintainer='Yanxian Li'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}

RUN apt-get update -y\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev\
&& apt-get install -y curl libcurl4-openssl-dev\
&& apt-get install -y git\
&& R -e "install.packages(c('BiocManager', 'remotes'))"\
&& R -e "BiocManager::install(c('dada2', 'biomformat', 'phyloseq', 'microbiome', 'philr', 'Maaslin2'))"\
&& R -e "remotes::install_github(c('rstudio/gt@efaee2d', 'jbisanz/MicrobeR@d542679','jbisanz/qiime2R@0.99.11'),upgrade='never')"

USER ${NB_USER}
RUN wget https://github.com/yanxianl/Li_AqFl2-Microbiota_ASM_2020/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2020-03-18/')); devtools::install_deps()"