FROM ubuntu:14.04
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN apt-get update -qq
RUN apt-get install -y software-properties-common
RUN apt-add-repository -y ppa:marutter/rrutter
RUN apt-add-repository -y ppa:marutter/c2d4u

RUN apt-get update -qq
RUN apt-get install -y --no-install-recommends r-base
RUN apt-get install -y r-cran-rcolorbrewer r-cran-ggplot2 
RUN apt-get install -y r-cran-igraph 
RUN apt-get install -y r-cran-vegan r-cran-fields
RUN apt-get install -y r-cran-pander

ADD root/PEDvweeklyreport-state-ts-01-08-14.csv /root/PEDvweeklyreport-state-ts-01-08-14.csv
ADD root/edge-proportionalities.csv /root/edge-proportionalities.csv
ADD root/mantel-testing.r /root/mantel-testing.r
ADD root/state_neighbors_fips.txt /root/state_neighbors_fips.txt 
ADD root/run-analysis /root/run-analysis

# Set environment variables.
ENV HOME /root

# Define working directory.
WORKDIR /root

CMD ["bash"]
