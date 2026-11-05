# RNAdecayCafe_paper_code
Scripts to reproduce figures in RNAdecayCafe manuscript

Required data:

1. RNAdecayCafe_v1.1_onetable.csv from the [Zenodo repo](https://zenodo.org/records/15785218?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjUzZTk5N2Y2LTdiYjAtNDVkYi05ZWI1LTU3ODgzNjY0ODMwMiIsImRhdGEiOnt9LCJyYW5kb20iOiJjNTZlOGQwNDM3ZjI1NDViZDQxYjBmOGEwN2IxYTcwYSJ9._RcoAZUq3dmNklIorm5yw9w87RuWEV75P-9quWkqHjpaCK_T2mhnCjQs_3J1IjQQ4gIgfrzE_j5Cyhm_4LDrQg).
2. The table of half-life/kdeg estimates from [Agarwal et al. 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02811-x#Sec22). This is 13059_2022_2811_MOESM2_ESM.xlsx in the Supplementary Information of that paper.
3. Table of all NR-seq derived half-life estimates.

2 and 3 can be obtained from the associated [Figshare repo](https://figshare.com/s/871ed87b4a845910b5be)

Required packages:

- data.table
- dplyr
- ggplot2
- tidyr
- readr
- purrr
- corrplot
- stringr
- pheatmap
- RcolorBrewer
- ComplexHeatmap
- circlize
- ggbreak
- broom
- intervals
