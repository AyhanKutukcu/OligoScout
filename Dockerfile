# OligoScout v0.4.1 — web engine + MCP server in one image.
#
# The image carries the code and the compiled C++ tools. Large data stay OUTSIDE
# the image and are mounted at /data (same layout as OligoScout/data):
#   /data/genomes/GRCh38.p14.fa                  (optional; needed only for index builds)
#   /data/indexes/ppfm/...                       (PPFM shards, built by RUN_SETUP)
#   /data/annotations/{chrM.fa, grch38_extra_contigs.fa/.tsv, *.tsv, gencode...gtf.gz}
#   and OligoScout/results/benchmarks/grch38_ppfm_batch/manifest.sa8.tsv -> /data/manifest
#
# Build:   docker build -t oligoscout:0.4.1 .
# Run MCP: docker run --rm -p 8766:8766 -v /srv/oligoscout-data:/data \
#            -e OLIGOSCOUT_MCP_TOKEN=... oligoscout:0.4.1
# Run web: docker run --rm -p 8765:8765 -v /srv/oligoscout-data:/data oligoscout:0.4.1 web

FROM ubuntu:24.04 AS build
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake ninja-build ca-certificates && \
    rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY OligoScout/ OligoScout/
RUN cd OligoScout && cmake --preset release >/dev/null && \
    cmake --build --preset release --target oligoscout_panel oligoscout_contig_scan \
          test_single_primer_products test_terminal_rule_search -j"$(nproc)" && \
    ./build/release/test_single_primer_products && ./build/release/test_terminal_rule_search

FROM ubuntu:24.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-numpy ca-certificates tini && \
    rm -rf /var/lib/apt/lists/* && \
    useradd --system --create-home --home-dir /opt/oligoscout oligoscout
WORKDIR /opt/oligoscout
COPY --chown=oligoscout:oligoscout reference_py/ reference_py/
COPY --chown=oligoscout:oligoscout web/ web/
COPY --chown=oligoscout:oligoscout mcp_server/ mcp_server/
COPY --chown=oligoscout:oligoscout LICENSE VERSION CHANGELOG.md ./
COPY --from=build --chown=oligoscout:oligoscout /src/OligoScout/build/release/oligoscout_panel \
     /src/OligoScout/build/release/oligoscout_contig_scan OligoScout/build/release/
RUN mkdir -p OligoScout/results/benchmarks/grch38_ppfm_batch && \
    ln -s /data OligoScout/data && \
    ln -s /data/manifest/manifest.sa8.tsv OligoScout/results/benchmarks/grch38_ppfm_batch/manifest.sa8.tsv && \
    cd reference_py && python3 -B -m unittest discover -s tests -q

COPY --chmod=755 deploy/docker-entrypoint.sh /usr/local/bin/oligoscout
USER oligoscout
ENV OLIGOSCOUT_ROOT=/opt/oligoscout \
    OLIGOSCOUT_GNOMAD_API=1 \
    PYTHONDONTWRITEBYTECODE=1
EXPOSE 8765 8766
VOLUME ["/data"]
ENTRYPOINT ["/usr/bin/tini", "--", "/usr/local/bin/oligoscout"]
CMD ["mcp"]
