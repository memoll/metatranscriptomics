# metatranscriptomics


To process data I used the SAMSA2 metatranscriptome analysis pipeline.

Documentation based on SAMSA2 tutorial: https://github.com/transcript/samsa2.

The statitiscal analyses presented here were all performed using R.

The present results are based on microbial annotations against the SEED Subsystems
hierarchical database (scripts starting with "ss") and the NCBI’s RefSeq bacterial genomes (scripts starting with "fb") and eukaryotic
genomes (scripts starting with "fu") databases. 

#### Steps of data analysis:
- DADA2 analysis of different sequening runs individually
- Merge the final results of multiple runs and assign taxonomy
- Quality control, denoising and cleaning data
- Explanatory analysis 
- Differential expression analysis 
