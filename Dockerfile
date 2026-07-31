FROM pytorch/pytorch:1.13.1-cuda11.6-cudnn8-runtime

LABEL org.opencontainers.image.title="GOBoost for FunctionBench/LAFA" \
      org.opencontainers.image.source="https://github.com/SunRui7/goboost-lafa" \
      org.opencontainers.image.description="Containerized GOBoost protein function predictor"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    GOBOOST_ESM_MODEL=/app/weights/esm1b_t33_650M_UR50S.pt

WORKDIR /app

ARG UBUNTU_MIRROR=archive.ubuntu.com
ARG UBUNTU_SECURITY_MIRROR=security.ubuntu.com
RUN sed -i \
      -e "s|archive.ubuntu.com|$UBUNTU_MIRROR|g" \
      -e "s|security.ubuntu.com|$UBUNTU_SECURITY_MIRROR|g" \
      /etc/apt/sources.list \
    && apt-get -o Acquire::Retries=5 update \
    && apt-get install -y --no-install-recommends ca-certificates unzip wget \
    && rm -rf /var/lib/apt/lists/*

COPY requirements-lafa.txt .
ARG PIP_INDEX_URL=https://pypi.org/simple
ARG PIP_TRUSTED_HOST=
ARG DGL_WHEEL_URL=https://data.dgl.ai/wheels/cu116/dgl-1.1.0%2Bcu116-cp310-cp310-manylinux1_x86_64.whl
RUN pip install --no-cache-dir -r requirements-lafa.txt \
    && wget --tries=5 --timeout=60 -q \
       -O '/tmp/dgl-1.1.0+cu116-cp310-cp310-manylinux1_x86_64.whl' "$DGL_WHEEL_URL" \
    && echo "00b222d27e9982356cab0a7c5a5101a728d95772ae2be1596ba4b7264a4fb3fa  /tmp/dgl-1.1.0+cu116-cp310-cp310-manylinux1_x86_64.whl" \
       | sha256sum -c - \
    && pip install --no-cache-dir \
       '/tmp/dgl-1.1.0+cu116-cp310-cp310-manylinux1_x86_64.whl' \
    && rm '/tmp/dgl-1.1.0+cu116-cp310-cp310-manylinux1_x86_64.whl'

# Freeze the GOBoost checkpoints in the image. The 7.29 GiB ESM-1b checkpoint
# is mounted separately, following LAFA's guidance for data larger than 5 GB.
ARG GOBOOST_MODEL_URL=https://zenodo.org/api/records/14048928/files/Model.zip/content
RUN mkdir -p /app/Model \
    && wget --tries=5 --timeout=60 --progress=dot:giga \
       -O /tmp/Model.zip "$GOBOOST_MODEL_URL" \
    && echo "bbdf045c13db349acab8b7777230184d  /tmp/Model.zip" | md5sum -c - \
    && unzip -q /tmp/Model.zip \
       'Model/best_*_All.pt' 'Model/best_*_Head.pt' 'Model/best_*_Tail.pt' \
       -d /app \
    && rm /tmp/Model.zip

COPY lafa_main.py .
COPY megraph ./megraph
COPY data/nrPDB-GO_2019.06.18_annot.tsv ./data/nrPDB-GO_2019.06.18_annot.tsv
COPY data/DMETrain/distribution_bp_All.txt ./data/DMETrain/distribution_bp_All.txt
COPY data/DMETrain/distribution_mf_All.txt ./data/DMETrain/distribution_mf_All.txt
COPY data/DMETrain/distribution_cc_All.txt ./data/DMETrain/distribution_cc_All.txt

ENTRYPOINT ["python", "/app/lafa_main.py"]
