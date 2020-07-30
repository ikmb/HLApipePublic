#!/bin/bash

BIM=$1


awk '{print "chr"$1":"$4"-"$4}' $BIM |sed 's/chr23/chrX/g'|sed 's/chr24/chrY/g' | sed 's/chr25/chrXY/g' | sed 's/chr26/chrMT/g' > prelift.pos
