# Augmented minimap2 with the tools the Sarek routes run in this image,
# plus leviosam2 0.5.0 and bamfix for the hap1-to-GRCh38 conversion.
# Build: docker build -t augmented-mapper-hap1 .
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc g++ make cmake libc6-dev zlib1g-dev libbz2-dev liblzma-dev libhts-dev libsdsl-dev \
        samtools bcftools tabix python3 python3-pysam mawk procps wget ca-certificates \
    && rm -rf /var/lib/apt/lists/*
COPY . /src
RUN make -C /src minimap2 && cp /src/minimap2 /usr/local/bin/ \
    && gcc -O2 -o /usr/local/bin/bamfix /src/tools/legacy_alignment/bamfix.c -lhts \
    && rm -rf /src
# leviosam2 as in bioinformatics_public/images (tests removed from CMakeLists.txt)
RUN wget -qO /tmp/leviosam2.tar.gz https://github.com/milkschen/leviosam2/archive/refs/tags/v0.5.0.tar.gz \
    && tar -xzf /tmp/leviosam2.tar.gz -C /tmp && cd /tmp/leviosam2-0.5.0 \
    && sed -i '209,289d' CMakeLists.txt \
    && printf '#ifndef VERSION\n#define VERSION "0.5.0"\n#endif\n' > src/version.hpp \
    && mkdir build && cd build && cmake -D CMAKE_BUILD_TYPE=release .. && make && make install \
    && rm -rf /tmp/leviosam2*
ENTRYPOINT []
