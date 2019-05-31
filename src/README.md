## apps to be used with genesis
To compile a genesis app, you need to add it to a local copy of [genesis](https://github.com/lczech/genesis) in the `apps` subfolder, then build the library using `make`. The resulting app can be found under `genesis/bin/apps`

### partial-tree-taxassign.cpp
Performs taxonomic assignment in a tree, where only some taxa have known taxonomic paths associated with them (in our case the rest were environmental sequences). In principle, this is achieved by first propagating taxonomic labels from teh tips toward the root, then checking which branches were not labeled (hence they lead to unlabaled/environmental taxa) and propagating the nearest label downwards.
