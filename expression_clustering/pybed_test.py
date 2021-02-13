from pybedtools import BedTool

eCRE = BedTool("reference/eCRE_locs/o2e33.bed")
archetypes = BedTool("reference/archetypes/mm10.archetype_motifs.v1.0.chr16.bed")

eCRE.intersect(archetypes,'-wb').saveas("reference/intersected.bed")