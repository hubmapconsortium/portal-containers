# Using Conda because pyarrow did not install easily on python base images.
FROM continuumio/miniconda3

# Copy requirements fits, so changes in main.py won't clear all cache layers.
COPY requirements.txt .
RUN pip install  -r ./requirements.txt

# Keep main installs, even if there are minor version changes in dependencies.
COPY requirements-freeze.txt .
RUN pip install  -r ./requirements-freeze.txt

COPY . .

CMD [ "python", "main.py", "--input_dir", "/input", "--output_dir", "/output" ]
