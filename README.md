export CONDA_ENVS_PATH=./.conda_litho_polyploid/envs

export CONDA_PKGS_DIRS=./.conda_litho_polyploid/pkgs

mamba env create --file resources/conda_env.yaml

## Install snakemake

mamba install -n litho_polyploid_env -c conda-forge -c bioconda snakemake=8.20
