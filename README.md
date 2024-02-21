######This repository contains the code and results for a basic nextflow workflow to trim and assemble paired end reads

Add all required packages to the yml file and initiate creation of a conda environment
```
conda env create -f nxtflow.yml
```

Activate the environment 
```
conda activate nts_hmwrk

```

Run the script provided in the repository
```
nextflow run nf_t_s.nf --reads '.fastq/*_{R1,R2}.fastq.gz'

```
