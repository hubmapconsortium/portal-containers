#! /bin/bash
# Removes all test-output-actual directories in the current directory
# Useful for local development when you want to re-run docker tests without
# having to manually delete the output directories.

sudo rm -rf `find . -type d -name test-output-actual`
