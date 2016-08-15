#!/bin/bash
# tracking production settings
#../../../mesord/mesord -i 1 -I 50 -c 1 -C -1 -E -p -g ligand_with_nucleoid_target.xml -t 60 -q 0.010 um -K -x 0.1 -w 0.0037

# fast settings for tuning the rates
../../../mesord/mesord -i 1 -I 50 -c 1 -C -1 -E -p -g ligand_with_nucleoid_target.xml -t 60 -q 0.050 um -K -x 0.1 -w 0.01
