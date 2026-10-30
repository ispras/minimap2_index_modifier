# Augmented minimap2 with the tools the Sarek route runs in this image.
# Build from the repository root: docker build -f docker/Dockerfile -t augmented-mapper .
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc make libc6-dev zlib1g-dev libhts-dev samtools bcftools tabix \
        python3 python3-pysam mawk procps \
    && rm -rf /var/lib/apt/lists/*
COPY . /src
RUN make -C /src minimap2 && cp /src/minimap2 /usr/local/bin/ && rm -rf /src
ENTRYPOINT []
