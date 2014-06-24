FROM ubuntu:14.04
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN apt-get update -qq
RUN apt-get install -y --no-install-recommends r-base
