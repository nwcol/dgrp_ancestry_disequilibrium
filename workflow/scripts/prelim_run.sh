#!/bin/bash






ancestry_hmm -i drgp_ahmm.input -s drgp_ahmm.samples -a 2 0.2 0.8 \
    -p 0 0.2 1598 -p 1 0.8 1598 -g -ne 25000000

