# We just need to use --file to point at it, instead of assuming it is in context.

# Using Conda because pyarrow did not install easily on python base images.
FROM --platform=linux/amd64 continuumio/miniconda3:26.3.2-2

# For tiff packages
RUN apt-get --allow-releaseinfo-change update && \
    apt-get install -y --no-install-recommends gcc python3-dev libhdf5-dev pkg-config python3-numcodecs && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --upgrade pip setuptools

# asciitree throwing a build error
RUN pip install --no-cache-dir --no-build-isolation asciitree

COPY requirements-freeze.txt .
RUN pip install --no-cache-dir -r ./requirements-freeze.txt

# In development, you may want to pin a single dependency in requirements.txt,
# without throwing away the entire cache layer from requirements-freeze.txt.
# (But once it works, you should check in an updated freeze!)

COPY requirements.txt .
RUN pip install --no-cache-dir -r ./requirements.txt

COPY . .

RUN useradd --create-home --uid 1000 appuser && \
    mkdir -p /input /output && \
    chown appuser:appuser /input /output
USER appuser

CMD [ "python", "main.py", \
      "--input_dir", "/input", \
      "--output_dir", "/output" ]
