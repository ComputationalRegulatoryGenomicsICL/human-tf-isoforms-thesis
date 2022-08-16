#!/usr/bin/env python

import sys
import re

inputYamlFile = sys.argv[1]
outputYamlFile = sys.argv[2]

with open(inputYamlFile, "r") as inputYaml, \
     open(outputYamlFile, "w") as outputYaml:
    prevSpaceCount = -1
    for line in inputYaml:
        spaceMatch = re.search('^[ ]+', line)
        if (spaceMatch):
            start, end = spaceMatch.span()
            currSpaceCount = end - start
        else:
            currSpaceCount = 0
        if (currSpaceCount <= prevSpaceCount):
            outputYaml.write(' ' * (prevSpaceCount + 2) + 'status: specific\n')
        outputYaml.write(line)
        prevSpaceCount = currSpaceCount
    outputYaml.write(' ' * (prevSpaceCount + 2) + 'status: specific\n')
