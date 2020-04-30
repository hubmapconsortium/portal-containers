#!/usr/bin/env bash
set -o errexit
set -o pipefail

die() { set +v; echo "$*" 1>&2 ; exit 1; }

while getopts "o:i:" opt
do
   case "$opt" in
      i ) INPUT_DIR="$OPTARG" ;;
      o ) OUTPUT_DIR="$OPTARG" ;;
   esac
done

for FILE in $INPUT_DIR/*ome.tif*
do
  # Handle both ways tiffs are named
  BASE_FILE_NAME=$(basename $FILE .ome.tif)
  BASE_FILE_NAME=$(basename $BASE_FILE_NAME .ome.tiff)

  N5_FILE=$OUTPUT_DIR/$BASE_FILE_NAME.n5/

  /opt/bioformats2raw/bin/bioformats2raw $FILE $N5_FILE  --tile_width 512 --tile_height 512 \
    || die "TIFF-to-n5 failed."
  echo "Wrote n5 pyramid from /input/$FILE to $N5_FILE"

  /opt/raw2ometiff/bin/raw2ometiff $N5_FILE $OUTPUT_DIR/$BASE_FILE_NAME.ome.tif --compression=zlib \
    || die "n5-to-pyramid failed."
  echo "Wrote OMETIFF pyramid from $N5_FILE to $OUTPUT_DIR/$FILE"
  
done
