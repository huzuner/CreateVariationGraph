# CreateVariationGraph
This repository contains code to create a variation graph needed to run several viral genome assemblers.

**Input**:
- Contig file: fasta format (filtered contigs)
- Forward reads: fastq format
- Reverse reads: fastq format
- Output directory: string
- (optional) --vg: path to the VG executable

**Output**:
- Variation graph: mod_graph.gfa
- Node abundances: node_abundances.txt

**Required dependencies**:
- minimap2 (bioconda)
- seqwish (bioconda)
- graph-tool (conda)
- VG (download executable from https://github.com/vgteam/vg/releases/tag/v1.64.1)
