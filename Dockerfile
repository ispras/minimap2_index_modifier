FROM phusion/baseimage:jammy-1.0.1

ARG HTSLIB_VERSION=1.17

ENV TERM=xterm-256color \
    HTSLIB_URL=https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2

COPY . minimap2_index_modifier

# install deps and cleanup apt garbage
RUN set -eux; \
  apt-get update; \
  apt-get install -y \
    wget \
    make \
    autoconf \
    gcc \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    bzip2 \
  && rm -rf /var/lib/apt/lists/*;

#install htslib
RUN set -eux; \
mkdir temp; \
cd temp; \
\
wget ${HTSLIB_URL}; \
tar -xf htslib-${HTSLIB_VERSION}.tar.bz2; \
cd htslib-${HTSLIB_VERSION}; \
\
autoreconf -i; \
./configure; \
make; \
make install; \
cd ../../; \
rm -rf temp;

#install minimap2_index_modifier
RUN set -eux; \
cd minimap2_index_modifier; \
make; \
cp minimap2 /usr/local/bin; \
cd ../; \
rm -rf minimap2_index_modifier;

ENV LD_LIBRARY_PATH=/usr/local/lib
