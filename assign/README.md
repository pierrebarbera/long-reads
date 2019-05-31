## results of taxonomic assignment
First level folder names specify what tree was used, and second level (where applicable) what set of queries were used.

Examples:
- For [pruned_long/long_read_queries](pruned_long/long_read_queries/) the taxonomic assignment was based on phylogenetic placement, with the reference tree being inferred from a combination of the full length queries, which were subsequently pruned out. The full length queries were then placed on that tree using `EPA-ng`, the result of which was fed into `gappa analyze assign`.

- For [comprehensive_long](comprehensive_long/) the tree infered from combined reference and full length query sequences was supplied to [`partial-tree-taxassign`](../src/partial-tree-taxassign.cpp)

The `V4` in this context means usage of the query alignments masked to only include the V4 hypervariable region.