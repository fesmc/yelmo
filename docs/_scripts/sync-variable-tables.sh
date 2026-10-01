#!/bin/sh
# Quarto pre-render step (see _quarto.yml): copy the variable io tables that
# the model reads (input/yelmo-variables-*.md) into the docs, so the docs
# tables cannot drift from the code. The copies are git-ignored.
set -e
cd "$(dirname "$0")/.."
for f in ../input/yelmo-variables-*.md; do
    cp "$f" .
done
