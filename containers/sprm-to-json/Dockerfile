# Using Conda because pyarrow did not install easily on python base images.
FROM continuumio/miniconda3:4.12.0

# For tiff packages
RUN apt-get --allow-releaseinfo-change update &&\
      apt-get install -y gcc python3-dev
COPY requirements-freeze.txt .
RUN pip install -r ./requirements-freeze.txt

# In development, you may want to pin a single dependency in requirements.txt,
# without throwing away the entire cache layer from requirements-freeze.txt.
# (But once it works, you should check in an updated freeze!)

COPY requirements.txt .
RUN pip install -r ./requirements.txt

COPY . .

CMD [ "python", "main.py", \
      "--input_dir", "/input", \
      "--output_dir", "/output" ]
