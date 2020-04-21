ARG BUILD_IMAGE=gradle:6.2.1-jdk8

#
# Build phase: Use the gradle image for building.
#
FROM ${BUILD_IMAGE} as build
USER root
RUN apt-get update -qq
RUN export TZ=America/New_York && export DEBIAN_FRONTEND=noninteractive && apt-get install -y tzdata
RUN apt-get install -y -qq zeroc-ice-all-runtime libblosc1
RUN mkdir /bioformats_pyramid && chown 1000:1000 /bioformats_pyramid
USER 1000

# Clone tiff-to-raw pyramid code.
WORKDIR /bioformats_pyramid
RUN git clone https://github.com/glencoesoftware/bioformats2raw.git

# Checkout and build.
WORKDIR /bioformats_pyramid/bioformats2raw
RUN git checkout tags/v0.1.0
RUN gradle build
RUN cd build/distributions && rm bioformats2raw*tar && unzip bioformats2raw*zip && rm -rf bioformats2raw*zip
USER root
RUN mv build/distributions/bioformats2raw* /opt/bioformats2raw
USER 1000
# Clone raw pyramid to tiff repo.
WORKDIR /bioformats_pyramid
RUN git clone https://github.com/glencoesoftware/raw2ometiff.git 

# Checkout and build.
WORKDIR /bioformats_pyramid/raw2ometiff
RUN git checkout tags/v0.1.0 
RUN gradle build
RUN cd build/distributions && rm raw2ometiff*tar && unzip raw2ometiff*zip && rm -rf raw2ometiff*zip
USER root
RUN mv build/distributions/raw2ometiff* /opt/raw2ometiff
USER 1000
# Set working directory containing new cli tools.
WORKDIR /opt
