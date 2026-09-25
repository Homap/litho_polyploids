[Repeated polyploidization shapes divergence in floral morphology in Lithophragma bolanderi (Saxifragaceae)] (https://www.pnas.org/doi/10.1073/pnas.2505119122)

**Abstract**
Polyploidization is an important driver of evolution and diversification in flowering plants. Here, we assess how repeated polyploidization may have shaped diversification of floral morphology in Lithophragma bolanderi (Saxifragaceae). This species comprises multiple cytotypes and varies geographically in its interactions with specialized pollinating moths in the genus Greya (Prodoxidae). Past studies have shown that coevolution with these moths has favored particular suites of floral characters but does not fully explain local and regional floral diversification. We combined phenotypic and genomic data from more than 1,800 individuals from 40 L. bolanderi populations spread across its entire range. Flow-cytometric analyses revealed a geographic mosaic of populations comprising one to four of three dominant (diploid, tetraploid, hexaploid) and three rare (triploid, pentaploid, octoploid) cytotypes. Whole-genome resequencing of a subset of populations suggested that polyploids arose from multiple autopolyploidization events, rather than a single event and/or through hybridization, albeit with some signals consistent with low levels of introgression from the congener Lithophragma glabrum. Quantification of flower traits from plants grown in a common garden showed that cytotype explained more than 15% of the variation in floral morphology, with polyploids showing more variability than diploids. Experimental induction of neopolyploids directly induced phenotypic changes but also indicated that local selection may have favored subsequent convergence in floral morphology among cytotypes in natural populations. Collectively, this comprehensive and integrative approach provides insights into how variability generating processes, such as polyploidization integrates with selection from species interactions to shape local floral diversification.

# Activate the environment

export CONDA_ENVS_PATH=./.conda_litho_polyploid/envs

export CONDA_PKGS_DIRS=./.conda_litho_polyploid/pkgs

mamba env create --file resources/conda_env.yaml

## Install snakemake

mamba install -n litho_polyploid_env -c conda-forge -c bioconda snakemake=8.20
