# Evaluation Strategy - Previous Evaluation


This evaluation strategy is performed in [1]. Li et al. do not perform CV on full data points, thereby avoiding the large number of negative instances during
evaluation. First, they create 10 folds of the data points belonging to group “1” for cross-validation. In each CV iteration, one fold from group “1” was set aside as the test set. Then, they randomly selected a subset of data points from group “0” equal in size to the fold from group “1”. Finally, two subsets of group “0” and group “1” are combined to generate the test set for that CV iteration. Thus, the majority of data points from group “0” are ignored in this evaluation because they are never included in the test sets. 

The idea is also demonstrated in ´´´crossval_method.m´´´ file, which can be found at https://github.com/lyhbio/HN-DRES/blob/main/Snakemake/BNNR/crossval_method.m


# Citation
[1] Li, Yinghong, Yinqi Yang, Zhuohao Tong, Yu Wang, Qin Mi, Mingze Bai, Guizhao Liang, Bo Li, and Kunxian Shu. "A comparative benchmarking and evaluation framework for heterogeneous network-based drug repositioning methods." Briefings in Bioinformatics 25, no. 3 (2024): bbae172.
