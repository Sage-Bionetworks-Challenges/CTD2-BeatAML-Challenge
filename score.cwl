#!/usr/bin/env cwl-runner
#
# Score SC1
#
cwlVersion: v1.0
class: CommandLineTool
baseCommand: score.py

hints:
  DockerRequirement:
    dockerPull: docker.synapse.org/syn20940521/scoring_harness:v1

inputs:
  - id: inputfile
    type: File
  - id: goldstandard
    type: File
  - id: check_validation_finished
    type: boolean?

arguments:
  - valueFrom: $(inputs.inputfile.path)
    prefix: -f
  - valueFrom: $(inputs.goldstandard.path)
    prefix: -g
  - valueFrom: results.json
    prefix: -r

requirements:
  - class: InlineJavascriptRequirement
     
outputs:
  - id: results
    type: File
    outputBinding:
      glob: results.json