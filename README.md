## Description
This tool modifies OrthoFinder's clustering pipeline by incorporating additional synteny constraints. During MCL graph construction, if two genes are linked by synteny (as defined in the .anchors file), an edge is automatically created between them in the similarity network, overriding the default BLAST score threshold requirement. All other steps remain consistent with standard OrthoFinder workflows.

## Usage
```
Orthofinder [-syn/--synteny <dir>] -f <dir> [option]
```
It accepts synteny files in the `.anchors` format (generated by *JCVI*).The file name should be like `SpeciesA.SpeciesB.anchors`
