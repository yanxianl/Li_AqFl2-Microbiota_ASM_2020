FROM rocker/binder:3.6.0
LABEL maintainer='Yanxian Li'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}

RUN apt-get update -y\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev\
&& apt-get install -y curl libcurl4-openssl-dev\
&& apt-get install -y git\
&& R -e "install.packages(c('BiocManager','devtools'))"\
&& R -e "BiocManager::install(c('biomformat','Biostrings','ComplexHeatmap','dada2','DECIPHER','microbiome','phangorn','phyloseq','philr'))"\
&& R -e "devtools::install_github(c('rstudio/gt@f793b33','jbisanz/MicrobeR@7207507','jbisanz/qiime2R@cd07f40','biobakery/Maaslin2@@RELEASE_3_10'),upgrade='never')"

USER ${NB_USER}

RUN wget https://github.com/yanxianl/Li_AqFl2-Microbiota_ASM_2020/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2020-03-16/')); devtools::install_deps()"

RUN rm DESCRIPTION.1; exit 0