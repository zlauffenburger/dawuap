FROM ubuntu:xenial

MAINTAINER dan.catalano@nwbt.co

RUN apt-get update && apt-get install -y python python-pip software-properties-common

RUN add-apt-repository ppa:ubuntugis/ppa; apt-get update; apt-get -y install nco gdal-bin libgdal-dev

COPY . /usr/src/app

WORKDIR /usr/src/app

RUN pip install -U pip pandas; python setup.py install

WORKDIR /usr/src/app/docs/example
