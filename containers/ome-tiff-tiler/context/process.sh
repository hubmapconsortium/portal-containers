#!/usr/bin/env bash

die() { set +v; echo "$*" 1>&2 ; exit 1; }

cd /input/

for FILE in *ome.tif*
do

  BASE_FILE_NAME=$(basename $FILE .ome.tif)
  N5_FILE=/output/$BASE_FILE_NAME.n5/

  SUCCESS=$(/opt/bioformats2raw/bin/bioformats2raw /input/$FILE $N5_FILE  --resolutions 8 --tile_width 512 --tile_height 512)
  if [[ $SUCCESS == 1 ]]; then
    die "TIFF-to-n5 failed."
  fi
  echo "Wrote n5 pyramid from /input/$FILE to $N5_FILE"

  SUCCESS=$(/opt/raw2ometiff/bin/raw2ometiff $N5_FILE /output/$FILE --compression=zlib)
  if [[ $SUCCESS == 1 ]]; then
    die "n5-to-pyramid failed."
  fi
  echo "Wrote OMETIFF pyramid from $N5_FILE to /output/$FILE"
  
done
