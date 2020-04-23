#!/usr/bin/env bash
cd /input/
for FILE in *ome.tif*
do
  BASE_FILE_NAME=$(basename $FILE .ome.tif)
  N5_FILE=/output/$BASE_FILE_NAME.n5/
  /opt/bioformats2raw/bin/bioformats2raw /input/$FILE $N5_FILE  --resolutions 8 --tile_width 512 --tile_height 512
  echo Wrote n5 pyramid from /input/$FILE to $N5_FILE
  /opt/raw2ometiff/bin/raw2ometiff $N5_FILE /output/$FILE --compression=zlib
  echo Wrote OMETIFF pyramid from $N5_FILE to /output/$FILE
done
