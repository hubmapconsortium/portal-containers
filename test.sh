#!/usr/bin/env bash
set -o errexit
set -o pipefail

./test_docker.sh
./test_cwl.sh
