Overview
======================

The trees

 - ref_and_long_reads_nws.fasta
 - ref_and_V4_reads_nws.fasta
 - ref_only_nws.fasta

were inferred with and without the read sequences.
Each dir also contains bootstrap calculations, 
as well as a "pruned" tree, where the read sequence tips
were removed again from the tree, so that they only
contain the ref taxa. That is, for the ref only tree,
no taxa were removed.

These pruned trees were then compared in `rf-dists`
to see where they differ from each other.
