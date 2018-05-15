FROM biobakery/picrust:latest
MAINTAINER Sam Minot sminot@fredhutch.org

# Small utility to facilitate execution via sciluigi
RUN pip install bucket_command_wrapper==0.3.0

# R dependencies
user root
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('ape')"
user linuxbrew
