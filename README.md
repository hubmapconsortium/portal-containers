# portal-containers

Docker containers to pre-process data for visualization in the portal.

The subdirectories in this repo all have the same structure:

- `context/`: A Docker context, including a `Dockerfile` and typically
  `main.py`, `requirements.txt`, and `requirements-freeze.txt`.
- `test-input/`, `test-output-actual/`, `test-output-expected/`: Test fixtures.
- `VERSION`: contains a semantic version number
- and a `README.md`.

Images are named by the containing directory.
Running `test.sh` will build (and test!) all the images.
You can then define `$INPUT_DIR`, `$OUTPUT_DIR`, and `$IMAGE`
to run an image with your own data:

```
docker run \
  --mount type=bind,source=$INPUT_DIR,target=/input \
  --mount type=bind,source=$OUTPUT_DIR,target=/output \
  $IMAGE
```

To push the latest versions to dockerhub just run:

```
test_docker.sh push
```

## Getting it to run in production

This repo is included as a submodule in [`ingest-pipeline`](https://github.com/hubmapconsortium/ingest-pipeline/tree/master/src/ingest-pipeline/airflow/dags/cwl):
When there are changes here that you want run in production, make a PR there to update that submodule to the latest code here, and make Joel a reviewer on the PR. Depending on the rate of change, it might be good to have a weekly routine of making PRs to `ingest-pipeline`. TBD.

In addition, each workflow must have a corresponding `-manifest.json` file conforming to [this schema](https://github.com/hubmapconsortium/ingest-pipeline/blob/devel/src/ingest-pipeline/schemata/pipeline_file_manifest.json), which has a `pattern`, `description`, and `edam_ontology_term` entry for each output file. (see [here](http://edamontology.org/page) for information about EDAM).

```
# In ingest-pipeline:
git checkout devel
git pull
git submodule update --init --recursive # This fails right now because we're not using plain URLs in .gitmodules.
git checkout -b username/update-portal-containers
cd src/ingest-pipeline/airflow/dags/cwl/portal-containers/
git checkout master
git pull
cd -
git commit -am 'Update portal-containers'
git push origin
# And then make a PR at: https://github.com/hubmapconsortium/ingest-pipeline
```

Here is a template for the PR into `ingest-pipelines`, if there is a new pipeline in `portal-containers` that needs to be run.  This helps us communicate what the pipeline's input, output, and purpose are:

```
# --NAME OF THE PIPELINE--

## Input Pipeline/Original Dataset:

## Output Pipeline (Optional):

## Description:

```

For example:

```
# [ome-tiff-tiler](https://github.com/hubmapconsortium/portal-containers/blob/master/ome-tiff-tiler.cwl)

## Input Pipeline/Original Dataset + Files:

- High resolution imaging from Vanderbilt data (OME-TIFF) files, such as those in `/hive/hubmap/lz/Vanderbilt TMC/e6e9bb7c01d3cb9cdb31a8da857f8832/processedMicroscopy/`

## Output Pipeline:

- [ome-tiff-offsets](https://github.com/hubmapconsortium/portal-containers/blob/master/ome-tiff-offsets.cwl)

## Description:

This pipeline takes as input Vanderbilt's processed microscopy data and outputs an image pyramid for visualization.  In addition, the `ome-tiff-offsets` pipeline needs to be run the output of `ome-tiff-tiler` so that images with z-stacks/large numbers of channels can be efficiently visualized.

```