Comprehensive multi-omics analysis of tandem duplicator phenotypes in non-small cell lung cancer

The Tandem Duplicator Phenotype (TDP) stands out as a significant genomic feature commonly observed in non-small cell lung cancer (NSCLC). Comprehensive multi-omics datasets, encompassing DNA copy number variations, transcriptomics, somatic single-nucleotide variations, clinical records, and cell-line drug sensitivity profiles from projects like TCGA and CCLE, are meticulously examined to unveil the distinctive attributes associated with TDP.

# datasets

76 datasets was deposited in https://github.com/jacklee2thu/Comprehensive-multi-omics-analysis-of-tandem-duplicator-phenotypes-in-non-small-cell-lung-cancer/tree/main/data.
1. 300 TSGs and 250 OGs classified by Davoli et al.[Cumulative haploinsufficiency and triplosensitivity drive aneuploidy patterns and shape the cancer genome. Cell 155, 948-962, doi:10.1016/j.cell.2013.10.011 (2013)] for an in-depth analysis are deposited in TSG.txt and oncogene.txt;
2. CNV associtated files include Non_TDP_cnv.txt, group_1_cnv.txt, TCGA-LUAD.masked_cnv.txt, only_1_2_cnv.txt, only_1_3_cnv.txt, only_2_3_cnv.txt, only_2_cnv.txt, only_3_cnv.txt in LUAD. Non_TDP_cnv.txt, TCGA-LUSC.masked_cnv.txt, group_1_cnv.txt, only_1_2_cnv.txt, only_1_3_cnv.txt, only_1_cnv.txt, only_2_3_cnv.txt, only_2_cnv.txt, only_3_cnv.txt in LUSC;
3. somatic single-nucleotide variations associated files include only_2.mutations.txt, only_2_3.mutations.txt in LUAD, group_1.mutations.txt, only_2_3.mutations.txt in LUSC;
4. clinical records associated files include TCGA-LUAD.GDC_phenotype.txt, TCGA-LUSC.GDC_phenotype.txt;
5. transcriptomics associated files include https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Lung%20Squamous%20Cell%20Carcinoma%20(LUSC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443;
6. cell-line drug sensitivity profiles associated files include LUAD_sample.txt, LUNG_drug.txt, LUSC_sample.txt, lung_sample.txt, and https://xenabrowser.net/datapages/?cohort=Cancer%20Cell%20Line%20Encyclopedia%20(CCLE)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443.

# script parameters

The scirpts was deposited in https://github.com/jacklee2thu/Comprehensive-multi-omics-analysis-of-tandem-duplicator-phenotypes-in-non-small-cell-lung-cancer/tree/main/script
The scripts LUSC.R and LUAD.R include 5 parameters, working_dictory indicates your working path; span_size_length indicates the span size length of TD, the default is 1000;breaks_seq indicates the interval size of TD, the default is seq(-1,6,0.1);fold_chang indicates the fold change cutoff of upregulated or downregulated genes; p_value indicates the p value cutoff of upregulated or downregulated genes (adj_p indicates the adjust p value cutoff of upregulated or downregulated genes in LUAD.R).
The scripts LUSC.R and LUAD.R also include the outcome files figure 1-6.

The outcome files are deposited in result directory.
