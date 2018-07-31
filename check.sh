#!/bin/bash

echo "TREES: "
var="$(find trees -name *.bestTree | tee /dev/tty | wc -l)"
echo "$var / 6"
echo "BOOTSTRAPS: "
var="$(find trees -name *.bootstraps | tee /dev/tty | wc -l)"
echo "$var / 30"

