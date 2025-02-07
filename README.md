# Germline pathogenic variation impacts somatic alterations and patient outcomes in pediatric CNS tumors

Ryan J. Corbett^, Rebecca Kaufman^, Shelly W. McQuaid, Zalman Vaksman, Saksham Phul, Miguel A. Brown, Jennifer L. Mason, Sebastian M. Waszak, Bo Zhang, Chuwei Zhong, Heena Desai, Ryan Hausler, Ammar S. Naqvi, Antonia Chroni, Zhuangzhuang Geng, Elizabeth M. Gonzalez, Yuankun Zhu, Allison P. Heath, Marilyn Li, Penn Medicine BioBank, Regeneron Genetics Center, Phillip B. Storm, Adam C. Resnick, Kara N. Maxwell, Kristina A. Cole, Angela J. Waanders, Miriam Bornhorst, Suzanne MacFarland, Jo Lynne Rokita+, Sharon J. Diskin+

^Equal authorship

+Co-senior authors

:tada: This work is now preprinted on [MexRxiv](https://www.medrxiv.org/content/10.1101/2025.02.04.25321499v1). :newspaper:


__DISCLAIMER__: This repository contains de-identified data only, and does NOT contain any raw germline variant data.  

## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct analyses for the manuscript noted above.

1. Clone the repository
```
git clone git@github.com:diskin-lab-chop/pbta-germline-somatic.git
```

2. Pull the docker container:
```
docker pull pgc-images.sbgenomics.com/diskin-lab/pbta-germline:1.0.0
```

3. Start the docker container

From the `pbta-germline-somatic` folder, run:

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/pbta-germline-somatic pgc-images.sbgenomics.com/diskin-lab/pbta-germline:1.0.0
```
 
Users can also run Rstudio in the project docker container from a web browser using the instructions below:

__Local Development in Rstudio__ (Max OS X and Linux users only)

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/pbta-germline-somatic pgc-images.sbgenomics.com/diskin-lab/pbta-germline:1.0.0
```

Then, navigate to `localhost:8787` in your web browser. The username for login is `rstudio` and the password will be whatever password is set in the `docker run` command above (default: `pass`).

__Development using Amazon EC2__

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 80:8787 -v $PWD:/home/rstudio/pbta-germline-somatic pgc-images.sbgenomics.com/diskin-lab/pbta-germline:1.0.0
```

Then, paste the instance IP address into your browser to start Rstudio. 

4. To execute shell within the docker image, from the `pbta-germline-somatic` folder, run:
```
docker exec -ti <CONTAINER_NAME> bash
```

5. Run the `download-data.sh` shell script to obtain latest data files: 
```
bash download_data.sh
```

6. Navigate to an analysis module and run the shell script:
```
cd /home/rstudio/pbta-germline-somatic/analyses/module_of_interest
```


### Below is the main directory structure listing the analyses and data files used in this repository

```
.
├── Dockerfile
├── README.md
├── Rplots.pdf
├── analyses
│   ├── alternative-splicing
│   ├── bed-intersect
│   ├── collapse-tumor-histologies
│   ├── cpg-enrichment
│   ├── demo-clin-stats
│   ├── dna-repair-variant-summary
│   ├── gene-expression
│   ├── germline-sv
│   ├── hgnc-liftover
│   ├── methylation
│   ├── mutational_signatures
│   ├── oncokb-annotation
│   ├── oncoprint
│   ├── predisposition-variants
│   ├── progression-analysis
│   ├── sample-distribution
│   ├── survival
│   ├── two-hits
│   └── variant-distribution
├── data
│   ├── cnv-cnvkit-cns.tsv.gz -> v9/cnv-cnvkit-cns.tsv.gz
│   ├── consensus_seg_annotated_cn_autosomes.tsv.gz -> v9/consensus_seg_annotated_cn_autosomes.tsv.gz
│   ├── consensus_seg_annotated_cn_x_and_y.tsv.gz -> v9/consensus_seg_annotated_cn_x_and_y.tsv.gz
│   ├── consensus_seg_with_status.tsv -> v9/consensus_seg_with_status.tsv
│   ├── cptac-protein-imputed-prot-expression-abundance.tsv.gz -> v9/cptac-protein-imputed-prot-expression-abundance.tsv.gz
│   ├── gene-counts-rsem-expected_count-collapsed.rds -> v9/gene-counts-rsem-expected_count-collapsed.rds
│   ├── gene-expression-rsem-tpm-collapsed.rds -> v9/gene-expression-rsem-tpm-collapsed.rds
│   ├── histologies-base.tsv -> v9/histologies-base.tsv
│   ├── histologies.tsv -> v9/histologies.tsv
│   ├── hope-protein-imputed-phospho-expression-abundance.tsv.gz -> v9/hope-protein-imputed-phospho-expression-abundance.tsv.gz
│   ├── hope-protein-imputed-prot-expression-abundance.tsv.gz -> v9/hope-protein-imputed-prot-expression-abundance.tsv.gz
│   ├── independent-specimens.methyl.primary-plus.tsv -> v9/independent-specimens.methyl.primary-plus.tsv
│   ├── independent-specimens.methyl.primary.tsv -> v9/independent-specimens.methyl.primary.tsv
│   ├── independent-specimens.rnaseqpanel.primary-plus.tsv -> v9/independent-specimens.rnaseqpanel.primary-plus.tsv
│   ├── independent-specimens.rnaseqpanel.primary.tsv -> v9/independent-specimens.rnaseqpanel.primary.tsv
│   ├── independent-specimens.wgs.primary-plus.tsv -> v9/independent-specimens.wgs.primary-plus.tsv
│   ├── independent-specimens.wgs.primary.tsv -> v9/independent-specimens.wgs.primary.tsv
│   ├── independent-specimens.wgswxs.primary-plus.tsv -> v9/independent-specimens.wgswxs.primary-plus.tsv
│   ├── independent-specimens.wgswxs.primary.tsv -> v9/independent-specimens.wgswxs.primary.tsv
│   ├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv -> v9/independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
│   ├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv -> v9/independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
│   ├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv -> v9/independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
│   ├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv -> v9/independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
│   ├── pbta-all-gene-loh.tsv.gz -> v9/pbta-all-gene-loh.tsv.gz
│   ├── pbta-cnv-consensus-gistic.zip -> v9/pbta-cnv-consensus-gistic.zip
│   ├── pbta-cnv-consensus.seg.gz -> v9/pbta-cnv-consensus.seg.gz
│   ├── pbta-fusion-putative-oncogenic.tsv -> v9/pbta-fusion-putative-oncogenic.tsv
│   ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds -> v9/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
│   ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds -> v9/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
│   ├── pbta-merged-plp-variants-autogvp-abridged-lowVAF-cpgs.tsv -> v9/pbta-merged-plp-variants-autogvp-abridged-lowVAF-cpgs.tsv
│   ├── pbta-merged-plp-variants-autogvp-abridged.tsv -> v9/pbta-merged-plp-variants-autogvp-abridged.tsv
│   ├── pbta-merged-plp-variants-autogvp-full.tsv -> v9/pbta-merged-plp-variants-autogvp-full.tsv
│   ├── pbta-plp-all-loh.tsv -> v9/pbta-plp-all-loh.tsv
│   ├── pbta-sv-manta.tsv.gz -> v9/pbta-sv-manta.tsv.gz
│   ├── pbta_germline_svs.tsv -> v9/pbta_germline_svs.tsv
│   ├── release-notes.md -> v9/release-notes.md
│   ├── snv-consensus-plus-hotspots.maf.tsv.gz -> v9/snv-consensus-plus-hotspots.maf.tsv.gz
│   ├── snv-mutation-tmb-coding.tsv -> v9/snv-mutation-tmb-coding.tsv
│   ├── splice-events-rmats.tsv.gz -> v9/splice-events-rmats.tsv.gz
│   ├── v2
│   ├── v3
│   ├── v4
│   ├── v5
│   ├── v6
│   ├── v7
│   ├── v8
│   └── v9
├── doc
│   └── release-notes.md
├── download_data.sh
├── figures
│   └── theme.R
├── make_summary_figures.Rmd
├── make_summary_figures.nb.html
└── scripts
    ├── download-methyl.sh
    ├── install_bioc.r
    ├── install_github.r
    ├── run_analysis_modules.sh
    └── update_fusion_gene_symbols.py
```

## Code Authors

Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)) and Ryan Corbett ([@rjcorb](https://github.com/rjcorb))

## Contact

For questions, please submit an issue.

