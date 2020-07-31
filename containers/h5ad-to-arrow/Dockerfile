# We just need to use --file to point at it, instead of assuming it is in context.

# Using Conda because pyarrow did not install easily on python base images.
FROM continuumio/miniconda3:4.7.12

# For tiff packages
RUN apt-get update &&\
      apt-get install -y gcc python3-dev
COPY requirements-freeze.txt .
RUN pip install  -r ./requirements-freeze.txt

# Download the GRCh38 Ensembl release to convert gene IDs to gene names.
RUN mkdir -p /cache/pyensembl
ENV PYENSEMBL_CACHE_DIR /cache/pyensembl
ENV ENSEMBL_RELEASE 99
RUN pyensembl install \
      --release $ENSEMBL_RELEASE \
      --species human \
      --custom-mirror http://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/

# In development, you may want to pin a single dependency in requirements.txt,
# without throwing away the entire cache layer from requirements-freeze.txt.
# (But once it works, you should check in an updated freeze!)

COPY requirements.txt .
RUN pip install  -r ./requirements.txt

COPY . .

CMD [ "python", "main.py", \
      "--input_dir", "/input", \
      "--output_dir", "/output" ]
