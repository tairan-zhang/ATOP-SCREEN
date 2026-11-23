FROM rocker/r-ver:4.6.1@sha256:780bc87219250a6dfffef6c812999731986fd94d8d5f7e572d1157c923d1b4bb
# The Linux image supplies the analysis worker and both compiled engines.
WORKDIR /opt/atop
# Resolve system packages against the same Ubuntu archive on every build.
RUN . /etc/os-release && test "$ID" = ubuntu && \
    rm -f /etc/apt/sources.list.d/*.list /etc/apt/sources.list.d/*.sources && \
    printf 'deb [check-valid-until=no] https://snapshot.ubuntu.com/ubuntu/20261101T000000Z/ %s main universe\n' \
      "$VERSION_CODENAME" "$VERSION_CODENAME-updates" "$VERSION_CODENAME-security" >/etc/apt/sources.list && \
    apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential gfortran cmake pkg-config ca-certificates git curl \
    libcurl4-openssl-dev libuv1-dev libnlopt-dev libssl-dev libxml2-dev libicu-dev libsqlite3-dev \
    libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
    libpng-dev libjpeg-dev libtiff-dev libwebp-dev libcairo2-dev libglpk-dev \
    libgmp-dev libgsl-dev libbz2-dev liblzma-dev zlib1g-dev libopenblas-dev \
    python3 python3-numpy python3-scipy && rm -rf /var/lib/apt/lists/*
ENV R_LIBS=/opt/atop/library ATOP_ANALYSIS_ROOT=/opt/atop/packages/analysis-engine \
    ATOP_CPP_ENGINE_LIBRARY=/opt/atop/lib/atop_permutation_engine.so ATOP_MAGECK_BIN=/opt/atop/mageck/bin
COPY renv.lock .
COPY renv/cellar/ renv/cellar/
COPY packaging/sources packaging/sources
COPY scripts/restore-r-runtime.R scripts/restore-r-runtime.R
ARG ATOP_BUILD_JOBS=2
RUN MAKEFLAGS=-j2 ATOP_BUILD_JOBS="$ATOP_BUILD_JOBS" Rscript --vanilla scripts/restore-r-runtime.R --library=/opt/atop/library
RUN mkdir -p /opt/atop/mageck/src /opt/atop/mageck/bin && \
    tar -xzf packaging/sources/mageck-0.5.9.5.tar.gz --strip-components=1 -C /opt/atop/mageck/src && \
    cd /opt/atop/mageck/src && make -C rra && make -C gsea && \
    cp bin/RRA bin/mageckGSEA /opt/atop/mageck/bin && \
    cp bin/mageck /opt/atop/mageck/bin/mageck && \
    sed -i '/^from __future__ import print_function/a import sys; sys.path.insert(0, "/opt/atop/mageck/src")' /opt/atop/mageck/bin/mageck && \
    ln -s /usr/bin/python3 /opt/atop/mageck/bin/python && \
    test "$(/opt/atop/mageck/bin/python -IB /opt/atop/mageck/bin/mageck --version)" = 0.5.9.5
COPY packages/analysis-engine packages/analysis-engine
COPY scripts/build-macos-cpp-engine.R scripts/build-macos-cpp-engine.R
RUN Rscript --vanilla scripts/build-macos-cpp-engine.R --output=/opt/atop/lib/atop_permutation_engine.so
COPY scripts/check-worker-designs.py scripts/check-worker-designs.py
RUN python3 scripts/check-worker-designs.py
ENTRYPOINT ["Rscript","--vanilla","/opt/atop/packages/analysis-engine/inst/worker/analysis_worker.R"]
