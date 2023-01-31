# controlFreq
The elaborate method for quantification of allele-specific expression. 

This R-package performs **sample overdispersion estimation** from RNA spike-ins and technical replication. 

## Installation:
```
devtools::install_github("gimelbrantlab/controlFreq")
```

## Example usage:
```
compute_iQCC_for_selected_samples(df = allelic-counts-table, reps = 1:10, sup_Q = 30)
```

* For allelic bias test and differential AI test, see [Qllelic R-package](https://github.com/gimelbrantlab/Qllelic) functions and *["Replicate sequencing libraries are important for quantification of allelic imbalance", A.Mendelevich et.al](https://www.nature.com/articles/s41467-021-23544-8)*.

* For data preprocessing procedure used in the paper and relevant recommendations, see [fastq2allelictabs](https://github.com/gimelbrantlab/fastq2allelictabs) repository.

![pic](https://github.com/gimelbrantlab/fastq2allelictabs/blob/main/schemes/ControlFreq_for_GitHub.png)

Figure was made in biorender.
