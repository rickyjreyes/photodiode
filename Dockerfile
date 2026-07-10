FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

COPY requirements.txt ./
RUN python -m pip install --upgrade pip \
    && python -m pip install --requirement requirements.txt

COPY . /app

RUN mkdir -p /work/outputs_npy_logcos
WORKDIR /work

VOLUME ["/work/outputs_npy_logcos"]

CMD ["python", "/app/scan_npy_logcos.py", "--input", "/app/captures/20260309_010802/waveform_raw.npy", "--sample-rate", "20000", "--fmin", "20", "--fmax", "1000", "--nperm", "100", "--cpu"]
