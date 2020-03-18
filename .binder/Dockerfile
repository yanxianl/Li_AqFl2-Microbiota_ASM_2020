FROM rocker/binder:3.6.0
LABEL maintainer='Yanxian Li'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}
RUN wget https://github.com/yanxianl/Li_AqFl2-Microbiota_ASM_2020/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2020-03-17/')); devtools::install_deps()"
