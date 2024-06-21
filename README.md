Comprehensive multi-omics analysis of tandem duplicator phenotypes in non-small cell lung cancer

The Tandem Duplicator Phenotype (TDP) stands out as a significant genomic feature commonly observed in non-small cell lung cancer (NSCLC). Comprehensive multi-omics datasets, encompassing DNA copy number variations, transcriptomics, somatic single-nucleotide variations, clinical records, and cell-line drug sensitivity profiles from projects like TCGA and CCLE, are meticulously examined to unveil the distinctive attributes associated with TDP.

69 datasets was deposited in https://github.com/jacklee2thu/Comprehensive-multi-omics-analysis-of-tandem-duplicator-phenotypes-in-non-small-cell-lung-cancer/tree/main/data.
The scirpts was deposited in https://github.com/jacklee2thu/Comprehensive-multi-omics-analysis-of-tandem-duplicator-phenotypes-in-non-small-cell-lung-cancer/tree/main/script

The input example data includes DNA copy number variations (TCGA-LUAD.masked_cnv.txt,TCGA-LUSC.masked_cnv.txt), transcriptomics (https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Lung%20Squamous%20Cell%20Carcinoma%20(LUSC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443), somatic single-nucleotide variations (LUSC_only_2_3.maf), clinical records (TCGA-LUAD.GDC_phenotype.txt, TCGA-LUSC.GDC_phenotype.txt), and cell-line drug sensitivity profiles (https://xenabrowser.net/datapages/?cohort=Cancer%20Cell%20Line%20Encyclopedia%20(CCLE)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).

The scripts LUSC.R and LUAD.R include 5 parameters, working_dictory indicates your working path; span_size_length indicates the span size length of TD, the default is 1000;breaks_seq indicates the interval size of TD, the default is seq(-1,6,0.1);fold_chang indicates the fold change cutoff of upregulated or downregulated genes; p_value indicates the p value cutoff of upregulated or downregulated genes (adj_p indicates the adjust p value cutoff of upregulated or downregulated genes in LUAD.R).
The scripts LUSC.R and LUAD.R also include the outcome files figure 1-6.

The outcome files are deposited in result directory.
