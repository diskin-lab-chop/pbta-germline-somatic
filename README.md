# Integrated germline and somatic profiling of pediatric CNS tumors reveals incidental findings and novel tumor biology

Ryan J. Corbett^, Rebecca Kaufman^, Shelly W. McQuaid, Zalman Vaksman, Saksham Phul, Miguel A. Brown, Jennifer L. Mason, Sebastian M. Waszak, Heena Desai, Ryan Hausler, Ammar S. Naqvi, Antonia Chroni, Zhuangzhuang Geng, Bo Zhang, Chuwei Zhong, Yuankun Zhu, Allison P. Heath, Marilyn Li, Phillip B. Storm, Adam C. Resnick, Kara N. Maxwell, Miriam Bornhorst, Kristina A. Cole, Angela J. Waanders, Suzanne MacFarland, Jo Lynne Rokita+, Sharon J. Diskin+

^Equal authorship

+Co-senior authors

## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct analyses for the manuscript noted above.

1. Clone the repository
```
git clone git@github.com:diskin-lab-chop/PBTA-germline.git
```

2. Pull the docker container:
```
docker pull pgc-images.sbgenomics.com/d3b-bixu/pbta-germline:latest
```

3. Start the docker container, from the `PBTA-germline` folder, run:
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/PBTA-germline pgc-images.sbgenomics.com/d3b-bixu/pbta-germline:latest
```

4. To execute shell within the docker image, from the `PBTA-germline` folder, run:
```
docker exec -ti <CONTAINER_NAME> bash
```

5. Run the `download-data.sh` shell script to obtain latest data files: 
```
bash download_data.sh
```

6. Navigate to an analysis module and run the shell script:
```
cd /home/rstudio/PBTA-germline/analyses/module_of_interest
```


### Below is the main directory structure listing the analyses and data files used in this repository

```
.
├── Dockerfile
├── README.md
├── analyses
│   ├── collapse-tumor-histologies
│   ├── cpg-enrichment
│   ├── demo-clin-stats
│   ├── dna-repair-variant-summary
│   ├── germline-sv
│   ├── hgnc-liftover
│   ├── oncokb-annotation
│   ├── predisposition-variants
│   ├── sample-distribution
│   ├── survival
│   ├── two-hits
│   └── variant-distribution
├── data
│   ├── consensus_seg_annotated_cn_autosomes.tsv.gz -> v6/consensus_seg_annotated_cn_autosomes.tsv.gz
│   ├── consensus_seg_annotated_cn_x_and_y.tsv.gz -> v6/consensus_seg_annotated_cn_x_and_y.tsv.gz
│   ├── consensus_seg_with_status.tsv -> v6/consensus_seg_with_status.tsv
│   ├── gene-expression-rsem-tpm-collapsed.rds -> v6/gene-expression-rsem-tpm-collapsed.rds
│   ├── histologies-base.tsv -> v6/histologies-base.tsv
│   ├── histologies.tsv -> v6/histologies.tsv
│   ├── independent-specimens.rnaseqpanel.primary-plus.tsv -> v6/independent-specimens.rnaseqpanel.primary-plus.tsv
│   ├── independent-specimens.rnaseqpanel.primary.tsv -> v6/independent-specimens.rnaseqpanel.primary.tsv
│   ├── independent-specimens.wgs.primary-plus.tsv -> v6/independent-specimens.wgs.primary-plus.tsv
│   ├── independent-specimens.wgs.primary.tsv -> v6/independent-specimens.wgs.primary.tsv
│   ├── independent-specimens.wgswxs.primary-plus.tsv -> v6/independent-specimens.wgswxs.primary-plus.tsv
│   ├── independent-specimens.wgswxs.primary.tsv -> v6/independent-specimens.wgswxs.primary.tsv
│   ├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv -> v6/independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
│   ├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv -> v6/independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
│   ├── pbta-all-gene-loh.tsv.gz -> v6/pbta-all-gene-loh.tsv.gz
│   ├── pbta-cnv-consensus-gistic.zip -> v6/pbta-cnv-consensus-gistic.zip
│   ├── pbta-cnv-consensus.seg.gz -> v6/pbta-cnv-consensus.seg.gz
│   ├── pbta-fusion-putative-oncogenic.tsv -> v6/pbta-fusion-putative-oncogenic.tsv
│   ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds -> v6/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
│   ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds -> v6/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
│   ├── pbta-plp-all-loh.tsv -> v6/pbta-plp-all-loh.tsv
│   ├── pbta-sv-manta.tsv.gz -> v6/pbta-sv-manta.tsv.gz
│   ├── pbta_germline_plp_calls_autogvp_abridged_11062023.tsv -> v6/pbta_germline_plp_calls_autogvp_abridged_11062023.tsv
│   ├── pbta_germline_plp_calls_autogvp_full_11062023.tsv -> v6/pbta_germline_plp_calls_autogvp_full_11062023.tsv
│   ├── release-notes.md -> v6/release-notes.md
│   ├── snv-consensus-plus-hotspots.maf.tsv.gz -> v6/snv-consensus-plus-hotspots.maf.tsv.gz
│   ├── snv-mutation-tmb-coding.tsv -> v6/snv-mutation-tmb-coding.tsv
│   └── v6
├── download_data.sh
├── figures
│   └── theme.R
└── scripts
    ├── install_bioc.r
    ├── install_github.r
    └── update_fusion_gene_symbols.py
```

## Code Authors

Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)) and Ryan Corbett ([@rjcorb](https://github.com/rjcorb))

## Contact

For questions, please submit an issue.

