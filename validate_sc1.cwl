#!/usr/bin/env cwl-runner
#
# Validate SC1
#
cwlVersion: v1.0
class: CommandLineTool
baseCommand: validate_sc1.py

hints:
  DockerRequirement:
    dockerPull: docker.synapse.org/syn20940521/scoring_harness:v1

inputs:
  - id: inputfile
    type: File
  - id: goldstandard
    type: File

arguments:
  - valueFrom: $(inputs.inputfile.path)
    prefix: -s
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

  - id: status
    type: string
    outputBinding:
      glob: results.json
      loadContents: true
      outputEval: $(JSON.parse(self[0].contents)['prediction_file_status'])

  - id: invalid_reasons
    type: string
    outputBinding:
      glob: results.json
      loadContents: true
      outputEval: $(JSON.parse(self[0].contents)['prediction_file_errors'])
