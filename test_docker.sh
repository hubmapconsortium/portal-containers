#!/usr/bin/env bash
set -o errexit
set -o pipefail

red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
reset=`tput sgr0`
start() { echo travis_fold':'start:$1; echo "$green$1$reset"; }
end() { set +v; echo travis_fold':'end:$1; echo; echo; }
die() { set +v; echo "$red$*$reset" 1>&2 ; exit 1; }

build_test() {
  TAG=$1
  BASENAME=$2
  docker build --file ./Dockerfile --tag $TAG context
  PWD_BASE=`basename $PWD`
  docker rm -f $PWD_BASE || echo "No container to stop"
  rm -rf test-output-actual || echo "No directory to delete"
  mkdir test-output-actual
  docker run \
    --name $PWD_BASE \
    --mount type=bind,source=$PWD/test-input/,target=/input \
    --mount type=bind,source=$PWD/test-output-actual/,target=/output \
    $TAG


  # hexdump -C test-output-expected/2x2.arrow > test-output-expected/2x2.arrow.hex.txt
  # hexdump -C test-output-actual/2x2.arrow > test-output-actual/2x2.arrow.hex.txt
  diff -w -r test-output-expected test-output-actual \
      --exclude=.DS_Store --exclude=*.arrow --exclude=*.ome.tif --exclude=*.ome.tiff \
      --exclude=*.ome.xml | head -n100 | cut -c 1-100

  # tiff-tiler is a special case: Java rather than Python.
  if [ "$BASENAME" != "ome-tiff-tiler" ]; then
    # `pip list --format=freeze` is the same as `pip freeze`
    # This works around a bug in pip freeze that creates a weird path instead of the package version:
    # https://stackoverflow.com/questions/62885911/pip-freeze-creates-some-weird-path-instead-of-the-package-version
    # This currently only occurs in the mudata-to-ui image since it uses a newer version of conda/pip.
    if [ "$BASENAME" == "mudata-to-ui" ]; then
      PIP_FREEZE="pip list --format=freeze"
    else
      PIP_FREEZE="pip freeze"
    fi
    diff <( docker run $TAG $PIP_FREEZE ) context/requirements-freeze.txt \
      || die "Update dependencies:
      docker run $TAG $PIP_FREEZE > $TAG/context/requirements-freeze.txt"
  fi

  MANIFEST_FILE=../../$BASENAME-manifest.json
  if [ ! -e "$MANIFEST_FILE" ]; then
    die "Manifest file does not exist for $MANIFEST_FILE"
  fi

  echo "$green$TAG is good!$reset"
}

for DIR in containers/*; do
  if [ -d "$DIR" ]; then
    pushd $DIR
      BASENAME=`basename $PWD`
      start $BASENAME
        VERSION=`cat VERSION`
        # Neither underscores nor double dash is allowed:
        # Don't get too creative!
        TAG="hubmap/portal-container-$BASENAME:$VERSION"
        build_test $TAG $BASENAME
        if [ "$1" == 'push' ]; then
          COMMAND="docker push $TAG"
          # If the current version has already been published, do not push.
          DOCKER_VERSIONS=`wget -q https://registry.hub.docker.com/v1/repositories/hubmap/portal-container-$BASENAME/tags -O - | jq -r '.[].name' \
            || die "No versions returned from dockerhub for $BASENAME.  If this is the intitial push for a container, please run \"$COMMAND\" before re-running this script."`
          CURRENT_VERSION=`cat VERSION`
          if grep -q "$CURRENT_VERSION" <<< "$DOCKER_VERSIONS"; then
            echo "$yellow Please update version of $BASENAME for pushing to docker.$reset"
          else
            echo "$green$COMMAND$reset"
            $COMMAND
          fi
        fi
      end $BASENAME
    popd
  fi
done
