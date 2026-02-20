#!/bin/bash
# . activate nextflow

#nextflow run /mnt/hdd2/anna/Myelo/src/project/annbeny/Myelo-Discovery -C nextflow.config -profile standard -resume -with-report
nextflow -C /mnt/hdd2/anna/Myelo/src/nextflow.config run /mnt/hdd2/anna/Myelo/src/project/annbeny/Myelo-Discovery -profile standard -resume -with-report