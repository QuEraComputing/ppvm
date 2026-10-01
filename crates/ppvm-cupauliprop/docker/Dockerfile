# Build environment for ppvm-cupauliprop: CUDA toolkit + cuPauliProp + Rust.
#
# Multi-arch (linux/amd64, linux/arm64). Without a GPU it builds, links and
# runs the host-only tests; with the NVIDIA container toolkit (`--gpus all`)
# it also runs the GPU tests and benchmarks. See check.sh.
FROM nvidia/cuda:12.8.1-devel-ubuntu24.04

ARG CUPAULIPROP_VERSION=0.5.0

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl unzip jq git build-essential \
    && rm -rf /var/lib/apt/lists/*

# cuPauliProp from the official PyPI wheel (just a zip with lib/ and include/).
RUN set -eux; \
    arch="$(uname -m)"; \
    url="$(curl -fsSL https://pypi.org/pypi/cupauliprop-cu12/${CUPAULIPROP_VERSION}/json \
        | jq -r --arg a "$arch" '.urls[] | select(.filename | endswith("manylinux2014_" + $a + ".whl")) | .url')"; \
    curl -fsSL -o /tmp/cupp.whl "$url"; \
    unzip -q /tmp/cupp.whl 'cuquantum/lib/*' 'cuquantum/include/*' -d /tmp/cupp; \
    mkdir -p /opt/cupauliprop; \
    mv /tmp/cupp/cuquantum/lib /tmp/cupp/cuquantum/include /opt/cupauliprop/; \
    rm -rf /tmp/cupp /tmp/cupp.whl

# Driver stub as libcuda.so.1 so binaries link and load on hosts without a
# GPU driver. Only on the path when NO_GPU=1 (see check.sh).
RUN mkdir -p /opt/cuda-stub \
    && ln -s /usr/local/cuda/lib64/stubs/libcuda.so /opt/cuda-stub/libcuda.so.1

ENV RUSTUP_HOME=/opt/rustup CARGO_HOME=/opt/cargo PATH=/opt/cargo/bin:$PATH
RUN curl -fsSL https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain none

ENV CUPAULIPROP_LIB_DIR=/opt/cupauliprop/lib \
    CUDA_HOME=/usr/local/cuda \
    LD_LIBRARY_PATH=/opt/cupauliprop/lib:/usr/local/cuda/lib64

WORKDIR /work
