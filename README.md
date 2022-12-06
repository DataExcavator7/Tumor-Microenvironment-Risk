![TMErisk](https://github.com/DataExcavator7/Tumor-Microenvironment-Risk/blob/Readme/TMErisk(1).png)
# Tumor-microenvironment Risk

Immune checkpoint inhibitor (ICI) treatment has brought the dawn of a new day for hepatocellular carcinoma (HCC) patients. However, only a minority of HCCs are benefitting from ICIs due to poor efficacy and safety concerns. As of today, few efficient predictive factors exist to precisely stratify the responders of immunotherapy in HCC. In this research, we developed a tumor microenvironment risk (TMErisk) model to divide HCC patients into different immune subtypes and considerably evaluate the prognosis of patients. Our result indicated that virally mediated HCC patients, that were more common in tumor protein P53 (TP53) alterations with lower TMErisk scores, were appropriate for ICI treatment, whereas multi-tyrosine kinase inhibitors (TKIs) were recommended for HCC patients with alcoholic hepatitis, that were more common in catenin beta 1 (CTNNB1) alterations with higher TMErisk scores.The TMErisk model in our research represents the first attempt to anticipate tumor tolerance in immune microenvironment via degree of immune infiltration in HCCs. We sincerely hope ICI treatment would be more precise and well-targeted for suitable populations in the near future.

The statistical analyses and visualization were performed by R software (version 4.0.5, except for oncoPredict package requesting R 4.2.0) and GraphPad Prism (version 8.0.1). All codes are storaged in [code](/example/profile.md).

We sincerely appreciate the researchers who provided algorithms used in this study, as well as carefully cite them in the text. The articles corresponding to the algorithm are listed as follows.

1.	Kolde, R., S. Laur, P. Adler, et al., Robust rank aggregation for gene list integration and meta-analysis. Bioinformatics, 2012. 28(4): p. 573-80.
2.	Yoshihara, K., M. Shahmoradgoli, E. Martínez, et al., Inferring tumour purity and stromal and immune cell admixture from expression data. Nat Commun, 2013. 4: p. 2612.
3.	Langfelder, P. and S. Horvath, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 2008. 9: p. 559.
4.	Langfelder, P. and S. Horvath, Fast R Functions for Robust Correlations and Hierarchical Clustering. J Stat Softw, 2012. 46(11).
5.	Bindea, G., B. Mlecnik, H. Hackl, et al., ClueGO: a Cytoscape plug-in to decipher functionally grouped gene ontology and pathway annotation networks. Bioinformatics, 2009. 25(8): p. 1091-3.
6.	Bindea, G., J. Galon, and B. Mlecnik, CluePedia Cytoscape plugin: pathway insights using integrated experimental and in silico data. Bioinformatics, 2013. 29(5): p. 661-3.
7.	Mayakonda, A., D.C. Lin, Y. Assenov, et al., Maftools: efficient and comprehensive analysis of somatic variants in cancer. Genome Res, 2018. 28(11): p. 1747-1756.
8.	Li, T., J. Fu, Z. Zeng, et al., TIMER2.0 for analysis of tumor-infiltrating immune cells. Nucleic Acids Res, 2020. 48(W1): p. W509-w514.
9.	Newman, A.M., C.L. Liu, M.R. Green, et al., Robust enumeration of cell subsets from tissue expression profiles. Nat Methods, 2015. 12(5): p. 453-7.
10.	Jiang, P., S. Gu, D. Pan, et al., Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response. Nat Med, 2018. 24(10): p. 1550-1558.
11.	Maeser, D., R.F. Gruener, and R.S. Huang, oncoPredict: an R package for predicting in vivo or cancer patient drug response and biomarkers from cell line screening data. Brief Bioinform, 2021. 22(6).
12.	Carvalho, B.S. and R.A. Irizarry, A framework for oligonucleotide microarray preprocessing. Bioinformatics, 2010. 26(19): p. 2363-7.
13.	Mehmood, A., A. Laiho, M.S. Venäläinen, et al., Systematic evaluation of differential splicing tools for RNA-seq studies. Brief Bioinform, 2020. 21(6): p. 2052-2065.
14.	Ritchie, M.E., B. Phipson, D. Wu, et al., limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res, 2015. 43(7): p. e47.
15.	Ginestet, C.J.J.o.t.R.S.S., ggplot2: Elegant Graphics for Data Analysis. 2011. 174(1): p. 245-246.
16.	Therneau, T.M. and P.M. Grambsch, Modeling Survival Data: Extending the Cox Model. 2013: Modeling Survival Data: Extending the Cox Model.
17.	Tibshirani, R., The lasso method for variable selection in the Cox model. Stat Med, 1997. 16(4): p. 385-95.
18.	Friedman, J., T. Hastie, and R. Tibshirani, Regularization Paths for Generalized Linear Models via Coordinate Descent. J Stat Softw, 2010. 33(1): p. 1-22.
19.	Blanche, P., J.F. Dartigues, and H. Jacqmin-Gadda, Estimating and comparing time-dependent areas under receiver operating characteristic curves for censored event times with competing risks. Stat Med, 2013. 32(30): p. 5381-97.
20.	Yu, G., L.G. Wang, Y. Han, et al., clusterProfiler: an R package for comparing biological themes among gene clusters. Omics, 2012. 16(5): p. 284-7.




