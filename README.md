# controlFreq
Pipeline for high-confidence allele-specific analysis of RNA-seq at scale. 

This R-package performs **sample overdispersion estimation** from RNA spike-ins or technical replication. 

Citation: A. Mendelevich, S. Gupta, A. Pakharev, A. Teodosiadis, A. A. Mironov, A. A. Gimelbrant. _Foreign RNA spike-ins enable accurate allele-specific expression analysis at scale_. Bioinformatics 39, i431–i439 (2023).
_[link](https://academic.oup.com/bioinformatics/article/39/Supplement_1/i431/7210455)_


## Installation:
```
docker pull sasha/ase
```
or
```
remotes::install_github("gimelbrantlab/controlFreq")
```

## Example usage:
```
compute_iQCC_for_selected_samples(df = allelic-counts-table, reps = 1:10, sup_Q = 30)
```
see [Wiki](https://github.com/gimelbrantlab/controlFreq/wiki) for more details.

## Related repositories:

* For allelic bias test and differential AI test, see [Qllelic R-package](https://github.com/gimelbrantlab/Qllelic) functions and *["Replicate sequencing libraries are important for quantification of allelic imbalance", A.Mendelevich et.al](https://www.nature.com/articles/s41467-021-23544-8)*.

* For data preprocessing procedure used in the paper and relevant recommendations, see [fastq2allelictabs](https://github.com/gimelbrantlab/fastq2allelictabs) repository.

![pic](https://github.com/gimelbrantlab/fastq2allelictabs/blob/main/schemes/ControlFreq_for_GitHub.png)

Figure above was made in [biorender](https://biorender.com/).
