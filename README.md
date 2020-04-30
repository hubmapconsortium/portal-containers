# portal-containers

Docker containers to pre-process data for visualization in the portal.

The subdirectories in this repo all have the same structure:

- `context/`: A Docker context, including at least
  `main.py`, `requirements.txt`, and `requirements-freeze.txt`.
- `test-input/`, `test-output-actual/`, `test-output-expected/`: Test fixtures.
- `VERSION`: contains a semantic version number
- and a `README.md`.

We have a single generic [`Dockerfile`](Dockerfile)
that we use to build an image from each context directory.
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
