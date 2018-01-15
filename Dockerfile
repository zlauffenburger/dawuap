FROM python:2.7-stretch

MAINTAINER dan.catalano@nwbt.co

RUN apt-get update && apt-get install -y nco && pip install pandas

COPY . /usr/src/app

WORKDIR /usr/src/app

RUN python setup.py install

WORKDIR /usr/src/app/docs/example
