# T2T-Polish

**Forked from main branch for modification & independent use** Evaluation and polishing workflows for T2T genome assemblies. 
Please cite if any of the codes shared in this repo was used:

Mc Cartney AM, Shafin K, Alonge M et al. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. Nat Methods (2022) doi: https://doi.org/10.1038/s41592-022-01440-3

For exact command lines and workflows used to generate the T2T-CHM13v1.0 and T2T-CHM13v1.1 assemblies, please refer to the [Methods](https://github.com/marbl/CHM13-issues#methods) section in the [CHM13-Issues](https://github.com/marbl/CHM13-issues) repo. Note that some of the tools have been updated since then, and are tracked on this repo.

### Automated polishing
* [Racon](https://github.com/isovic/racon/tree/liftover): Liftover branch for outputting edits in `.vcf`
* [Merfin](https://github.com/arangrhie/merfin): Latest stable code-base

### Base level QV estimation and K-Mer Databases
* [Meryl](https://github.com/marbl/meryl)
* [Merqury](https://github.com/marbl/merqury)
