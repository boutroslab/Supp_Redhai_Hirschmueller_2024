# get the base shiny app
FROM rocker/r-ver:4.3.1


# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libglpk-dev

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean
	
	
# create a directory where all important files will be located
RUN mkdir -p /home/decode_user/app
	
# set workdir
WORKDIR /home/decode_user/app
RUN cd /home/decode_user/app

# install renv
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# copy necessary files
COPY renv.lock renv.lock

# prepare package installation
RUN mkdir -p renv/library
ENV RENV_PATHS_LIBRARY renv/library


RUN mkdir renv/.cache
ENV RENV_PATHS_CACHE renv/.cache

# reinstall packages
RUN R -e "renv::restore()"

# copy necessary data
RUN mkdir data
COPY data/filtered_seurat_object.rds data/filtered_seurat_object.rds
COPY data/filtered_seurat_object_RNAi.rds data/filtered_seurat_object_RNAi.rds
COPY data/fbid_gene_mapping.csv data/fbid_gene_mapping.csv
COPY app.R app.R
COPY plot_theme.R plot_theme.R

# create a new user
RUN useradd -m decode_user -d /home/decode_user
USER decode_user

# expose port on which app will run. 
EXPOSE 5463

CMD ["R", "-e", "shiny::runApp('/home/decode_user/app/app.R', host = '0.0.0.0', port = 5463)"]
