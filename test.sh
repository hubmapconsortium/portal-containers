#!/usr/bin/env bash
set -o errexit
set -o pipefail

red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`
start() { echo travis_fold':'start:$1; echo "$green$1$reset"; }
end() { set +v; echo travis_fold':'end:$1; echo; echo; }
die() { set +v; echo "$red$*$reset" 1>&2 ; exit 1; }

test_docker() {
  NAME=`basename $PWD`
  start $NAME
  docker build --file ../Dockerfile --tag $NAME context
  docker rm -f ${NAME}-container || echo "No container to stop"
  rm -rf test-output-actual || echo "No directory to delete"
  mkdir test-output-actual
  docker run \
    --env TEXT_FOR_DIFF=true \
    --mount type=bind,source=$PWD/test-input/,target=/input \
    --mount type=bind,source=$PWD/test-output-actual/,target=/output \
    $NAME
  diff <( docker run $NAME pip freeze ) context/requirements-freeze.txt \
    || die "Update dependencies:
    docker run $NAME pip freeze > $NAME/context/requirements-freeze.txt"

  # hexdump -C test-output-expected/2x2.arrow > test-output-expected/2x2.arrow.hex.txt
  # hexdump -C test-output-actual/2x2.arrow > test-output-actual/2x2.arrow.hex.txt

  diff -w -r test-output-expected test-output-actual \
      --exclude=.DS_Store --exclude=*.arrow \
    | head -n100 | cut -c 1-100
  end $NAME
}

for DIR in *; do
  if [ -d "$DIR" ]; then
    pushd $DIR
      test_docker
    popd
  fi
done
