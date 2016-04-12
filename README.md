# snp-imputation-nf
Nextflow script to create the imputation

A portable, scalable snp imputation pipeline implemented in Nextflow.

### [Requirements](#requirements)

The pipeline is built on [Nextflow](http://nextflow.io) as a workflow engine, so it needs to be installed first:
```
curl -fsSL get.nextflow.io | bash
```

With Nextflow installed, the easiest way to use the pipeline is to use the prepared Docker container which contains all external dependencies.
```
docker pull insilicodb/docker-impute2
```

The imputation pipeline needs reference data this can be found on the [Impute2 website](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference). The reference data the pipeline is tested with is [here](https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz).

### [Running the pipeline](#running)

Here's how to start an example run using Docker (using the example dataset and parameterization included in the distribution):
```
$ nextflow run insilicodb/snp-imputation-nf -profile docker TODO ADD PARAMETERS!!!
```

For your own runs, provide your own file names, paths, parameters, etc. as defined in the `nextflow.config` file.

