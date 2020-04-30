# This is adapted in part from https://github.com/glencoesoftware/bioformats2raw/blob/v0.2.0/Dockerfile.

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

# Clone tiff-to-raw pyramid code.
WORKDIR /bioformats_pyramid
RUN git clone https://github.com/glencoesoftware/bioformats2raw.git

# Checkout and build.
WORKDIR /bioformats_pyramid/bioformats2raw
RUN git checkout tags/v0.1.1
RUN gradle build
RUN cd build/distributions && rm bioformats2raw*tar && unzip bioformats2raw*zip && rm -rf bioformats2raw*zip
RUN mv build/distributions/bioformats2raw* /opt/bioformats2raw

# Clone raw pyramid to tiff repo.
WORKDIR /bioformats_pyramid
RUN git clone https://github.com/glencoesoftware/raw2ometiff.git 

# Checkout and build.
WORKDIR /bioformats_pyramid/raw2ometiff
RUN git checkout tags/v0.1.3
RUN gradle build
RUN cd build/distributions && rm raw2ometiff*tar && unzip raw2ometiff*zip && rm -rf raw2ometiff*zip
RUN mv build/distributions/raw2ometiff* /opt/raw2ometiff
# Set working directory containing new cli tools.
WORKDIR /opt
COPY . /opt/
RUN chmod +x /opt/process.sh
CMD ["/opt/process.sh", "-o", "/output" , "-i", "/input"]