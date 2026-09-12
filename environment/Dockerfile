FROM ubuntu:24.04@sha256:33ceb71981b602c1a7443a53469e4dba065f7503eab3078a2d7a57a2ab987517
ARG DEBIAN_FRONTEND=noninteractive
COPY debs/r-base-core_4.5.3-1.2404.0_amd64.deb /tmp/r-base-core.deb
RUN echo '4066b325cc0e3e22ad7c786576fb8bef555ea0c09bd58d7a3f283969d87b9aa7  /tmp/r-base-core.deb' | sha256sum -c - \
 && apt-get update \
 && apt-get install -y --no-install-recommends /tmp/r-base-core.deb \
    build-essential gfortran libcurl4-openssl-dev libssl-dev libxml2-dev \
    python3 python3-venv ca-certificates \
 && rm -rf /var/lib/apt/lists/* /tmp/r-base-core.deb
COPY r-sources /tmp/r-sources
COPY r-packages.tsv r-source-lock.json install_packages.R verify_sources.py /tmp/
RUN python3 /tmp/verify_sources.py \
 && Rscript --vanilla /tmp/install_packages.R \
 && python3 -m venv /opt/repro-venv
COPY wheels /tmp/wheels
RUN /opt/repro-venv/bin/pip install --no-index --find-links=/tmp/wheels numpy==2.2.6 pandas==2.2.3 \
 && rm -rf /tmp/r-sources /tmp/wheels
ENV PATH=/opt/repro-venv/bin:$PATH
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0 TZ=UTC LC_ALL=C.UTF-8
WORKDIR /work
ENTRYPOINT []
CMD ["python3", "repro/run.py", "--help"]
