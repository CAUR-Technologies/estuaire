FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV SITE_DIR='/estuary/src/estuaire/site_scons'
ENV PATH=$PATH:'/estuary/src/eikonal-ng/bin'
ENV PYTHONPATH=$PYTHONPATH:'/estuary/src/agstd':$SITE_DIR

COPY . /estuary

RUN apt update \
    && apt update \
    && apt install python2 python2-dev make gcc-8 -y \
    && apt install vim -y \
    && apt install curl -y \
    && curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py \
    && apt remove curl -y \
    && python2 get-pip.py \
    && ln -s /usr/bin/python2.7 /usr/bin/python \
    && apt install build-essential autoconf libtool pkg-config -y \
    && pip install cython \
    && pip install mako numpy scipy \
    && pip install scons dateutils matplotlib tqdm pathlib \
    && cd /estuary/src/eikonal && make clean && make \
    && python setup.py install


