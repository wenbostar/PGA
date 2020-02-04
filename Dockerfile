FROM bioconductor/release_base2:R3.6.2_Bioc3.10


MAINTAINER wenbostar@gmail.com

RUN  rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*


RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get -y  install --fix-missing openjdk-8-jre \
    tcl8.6-dev \
    tk \
    expat \
    libexpat-dev && \
    apt-get clean





ADD install.R /tmp/

RUN R -f /tmp/install.R

RUN echo "R_LIBS=/usr/local/lib/R/host-site-library:\${R_LIBS}" > /usr/local/lib/R/etc/Renviron.site
RUN echo "R_LIBS_USER=''" >> /usr/local/lib/R/etc/Renviron.site
RUN echo "options(defaultPackages=c(getOption('defaultPackages'),'BiocManager'))" >> /usr/local/lib/R/etc/Rprofile.site

RUN chmod -R 755 /opt/

#change working directory
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]

