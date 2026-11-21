# Dockerfile for BBGT-GEM-Maxwell VLBI DSP Pipeline (v1.0-alpha)
FROM python:3.11-slim

LABEL maintainer="Matthias Asare-Darko <Darkon Research Laboratory, Ghana>"
LABEL description="Generalized Non-Invasive Quadrature Residual Filter & Planetary Kernel Demodulation Engine Anchored at GRAO (Kuntunse)"
LABEL version="1.0-alpha"
LABEL org.opencontainers.image.source="https://doi.org/10.5281/zenodo.24323557"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg

WORKDIR /opt/bbgt_vlbi

COPY requirements.txt /opt/bbgt_vlbi/
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

COPY gem_maxwell_pipeline.py /opt/bbgt_vlbi/
COPY uvfits_benchmark_and_diagnostics.py /opt/bbgt_vlbi/
COPY PROPOSAL.md /opt/bbgt_vlbi/
COPY README.md /opt/bbgt_vlbi/

RUN mkdir -p /opt/bbgt_vlbi/output
WORKDIR /opt/bbgt_vlbi/output

CMD ["sh", "-c", "python3 /opt/bbgt_vlbi/gem_maxwell_pipeline.py && python3 /opt/bbgt_vlbi/uvfits_benchmark_and_diagnostics.py"]
