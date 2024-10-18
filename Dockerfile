FROM rocker/shiny
LABEL maintainer="Vasileios Tsiamis <vasileios@bmb.sdu.dk>"
LABEL description="Docker image of ViQoR implementation on top of shiny-server. The number of to-be-installed R packages requires patience when building this image."

RUN apt-get update && apt-get install -y \
    libssl-dev \
    liblzma-dev \
    libbz2-dev \
    libicu-dev \
    libv8-dev \
    libglpk-dev \
    openjdk-8-jre && apt-get clean

#Packages
RUN R -e "install.packages('BiocManager', repos='http://cran.us.r-project.org'); \
          update.packages(ask=F)"
RUN R -e "library(BiocManager); \
          BiocManager::install(c('preprocessCore', 'Biostrings', 'shiny', 'shinydashboard', 'DT', 'protr', 'shinyjs', 'igraph', 'networkD3', \
                                 'plotly', 'shinyBS', 'shinycssloaders', 'ggplot2', 'reshape2', 'htmlwidgets', \
                                 'V8', 'seqinr', 'dplyr', 'seqinr', 'heatmaply', 'remotes'), ask=F)"
RUN R -e "library(remotes); \
          remotes::install_github('rstudio/webshot2')"

#Chrome
RUN wget "https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb"
RUN apt-get install -y ./google-chrome-stable_current_amd64.deb && rm google-chrome-stable_current_amd64.deb 
RUN sed -i '${s/$/ --no-sandbox/g}'  /opt/google/chrome/google-chrome

RUN rm -rf /srv/shiny-server
RUN mkdir /srv/shiny-server
RUN mkdir /srv/shiny-server/Dataset
RUN mkdir /srv/shiny-server/www
RUN mkdir /srv/shiny-server/smileys
RUN mkdir /srv/shiny-server/utilities-4.12.9
COPY *R  /srv/shiny-server/
COPY *pdf  /srv/shiny-server/
COPY Dataset/* /srv/shiny-server/Dataset/
COPY www/* /srv/shiny-server/www/
COPY smileys/* /srv/shiny-server/smileys/
COPY utilities-4.12.9/* /srv/shiny-server/utilities-4.12.9/



