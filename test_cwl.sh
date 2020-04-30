#!/usr/bin/env bash
set -o errexit
set -o pipefail

red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`
start() { echo travis_fold':'start:$1; echo "$green$1$reset"; }
end() { set +v; echo travis_fold':'end:$1; echo; echo; }
die() { set +v; echo "$red$*$reset" 1>&2 ; exit 1; }

OUTPUT_NAME=test-output-actual
for CWL_PATH in $PWD/*.cwl; do

  LABEL=`basename $CWL_PATH .cwl`
  CWL_FILE="$LABEL.cwl"
  start $LABEL

  cd $PWD/workflows/$LABEL/
  rm -rf $OUTPUT_NAME
  mkdir $OUTPUT_NAME
  # CWL needs be run inside the output folder.
  cd $OUTPUT_NAME
  ../../../$CWL_FILE ../test-job.yml
  # Go back to the root directory to look for other cwl files.
  cd ../../../

  diff -w -r ./workflows/$LABEL/test-output-expected ./workflows/$LABEL/$OUTPUT_NAME -x .DS_Store \
    -x *ome.xml.bak -x *.ome.tif -x *.ome.tiff -x *.ome.xml | head -n100 | cut -c 1-100

  end $LABEL
  
done