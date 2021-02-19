Nextflow pipelines for RNA-seq
==============================

What is Nextflow?
-----------------

[Nextflow](https://www.nextflow.io/) is system to run bioinformatics pipelines on HPC clusters such as _Deigo_.  Its main advantages are:

 - It can use batch job submission systems automatically (on _Deigo_, it is _SLURM_).
 - It downloads all the softwares it needs automatically.
 - It aggregates quality control information in a comprehensive report.
 - It enables reproducible and open research as the pipelines are archived for long-term and document precisely which versions of the softwares were used.
 
On _Deigo_, the `nextflow` software is provided as a module by the OIST Bioinformatics User Group ([BioinfoUgrp](https://github.com/oist/BioinfoUgrp)).
To load its most recent version, type `ml bioinfo-ugrp-modules Nextflow`.

_nf-core_ pipelies
------------------

The [nf-core](https://www.nextflow.io/) project provides high-quality pipelines for various kinds of quantitative sequencing data. 
In particular, it provides a [RNA-seq](https://nf-co.re/rnaseq) pipeline.

The _nf-core_ project also provides a tool of the same name (`nf-core`), to faciliate running the pipelines.
On _Deigo_, load its most recent version with `ml bioinfo-ugrp-modules nf-core`.

Running _nf-core_ pipelines on _Deigo_
--------------------------------------

Load the modules you need (see above) and go in your work directory on the `/flash` filesystem.

As Nextflow runs `sbatch` in the background, you can run the `nextflow` or the `nf-core` command from the login node.

After the computation is done, move its folder to the `/bucket` filesystem, as the space on `/flash` is limited.

By default, the final results are in the `results` subfolder.
The `work` subfolder contains some intermediary files, a copy of the final files, and a copy of the downloaded software
(as Singularity containers).  This folder may be useful to restart the computation in case it crashed or if you want
to change an option.  Otherwise it can be deleted.  If you would like to be extra cautious about reproducibility,
you can move the containers okt of the `work` folder before deleting it.
