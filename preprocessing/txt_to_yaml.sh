#!/usr/bin/env bash

<../../../data/interpro/ParentChildTreeFile.txt \
sed 's/:.*$/: /g' | \
sed 's/--/  /g' > \
../output/ParentChildTreeFile.yaml
