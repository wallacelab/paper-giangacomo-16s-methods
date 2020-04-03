#! /bin/bash

head -n 2 *.fasta --quiet > zymo_16s.fa  # .fa so it doesn't include itself if this gets rerun
