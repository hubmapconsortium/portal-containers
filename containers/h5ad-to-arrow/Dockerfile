# We just need to use --file to point at it, instead of assuming it is in context.

# Using Conda because pyarrow did not install easily on python base images.
FROM continuumio/miniconda3:4.7.12

# For tiff packages
RUN apt-get update &&\
      apt-get install -y gcc python3-dev
COPY requirements-freeze.txt .
RUN pip install  -r ./requirements-freeze.txt

# Download the GRCh38 Ensembl release 100 to convert gene IDs to gene names.
# Need to use custom URLs with http:// rather than ftp:// for travis,
# and becaues the pyensembl package only knows about releases up to 99.
RUN mkdir -p /cache
ENV PYENSEMBL_CACHE_DIR /cache
ENV ENSEMBL_URL http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
RUN pyensembl install \
      --species human \
      --reference-name GRCh38 \
      --annotation-name ensembl \
      --gtf $ENSEMBL_URL


# In development, you may want to pin a single dependency in requirements.txt,
# without throwing away the entire cache layer from requirements-freeze.txt.
# (But once it works, you should check in an updated freeze!)

COPY requirements.txt .
RUN pip install  -r ./requirements.txt

COPY . .

CMD [ "python", "main.py", \
      "--input_dir", "/input", \
      "--output_dir", "/output" ]
