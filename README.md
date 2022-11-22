# Germline variation across OpenPBTA patients

## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct analyses for the manuscript noted above.

1. Clone the repository
```
git clone https://github.com/diskin-lab-chop/OpenPBTA-germline.git
```

2. Pull the docker container:
```
docker pull pgc-images.sbgenomics.com/corbettr/pbta-germline:latest
```

3. Start the docker container, from the `OpenPBTA-germline` folder, run:
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/OpenPBTA-germline pgc-images.sbgenomics.com/corbettr/pbta-germline:latest
```

4. To execute shell within the docker image, from the `OpenPBTA-germline` folder, run:
```
docker exec -ti <CONTAINER_NAME> bash
```

5. Navigate to an analysis module and run the shell script:
```
cd /home/rstudio/OpenPBTA-germline/analyses/module_of_interest
```


### Below is the main directory structure listing the analyses and data files used in this repository

```
.
├── Dockerfile
├── README.md
├── analyses
│   ├── collapse-tumor-histologies
│   │   ├── 01-collapse-tumor-histologies.R
│   │   ├── input
│   │   │   └── PBTA_Germline_801_10102022.xlsx
│   │   └── results
│   │       ├── primary-tumor-histologies-collapsed-by-germline.tsv
│   │       ├── tumor-histologies-collapsed-by-germline-CBTN.tsv
│   │       └── tumor-histologies-collapsed-by-germline.tsv
│   ├── oncokb-annotation
│   │   ├── annotate-maf-oncokb.Rmd
│   │   ├── annotate-maf-oncokb.nb.html
│   │   ├── code
│   │   │   ├── 01-create-subset-maf.Rmd
│   │   │   ├── 01-create-subset-maf.nb.html
│   │   │   └── 02-run-oncokb-annotation.sh
│   │   ├── input
│   │   │   └── cpg.txt
│   │   ├── results
│   │   │   ├── snv-consensus-plus-hotspots-goi-oncokb.maf.tsv
│   │   │   └── snv-consensus-plus-hotspots-goi.maf.tsv
│   │   └── run-oncokb-annotation.sh
│   └── two-hits
│       ├── assess-two-hits.Rmd
│       ├── assess-two-hits.nb.html
│       └── results
│           ├── germline-somatic-collapsed-by-gene.tsv
│           ├── germline-somatic-two-gene-hits.tsv
│           └── pbta-oncokb-oncogenic-maf.tsv
├── data
│   └── v1
│       ├── consensus_seg_annotated_cn_autosomes.tsv.gz
│       ├── consensus_seg_annotated_cn_x_and_y.tsv.gz
│       ├── consensus_seg_with_status.tsv
│       ├── gene-expression-rsem-tpm-collapsed.rds
│       ├── independent-specimens.rnaseqpanel.primary-plus.tsv
│       ├── independent-specimens.rnaseqpanel.primary.tsv
│       ├── independent-specimens.wgs.primary-plus.tsv
│       ├── independent-specimens.wgs.primary.tsv
│       ├── independent-specimens.wgswxs.primary-plus.tsv
│       ├── independent-specimens.wgswxs.primary.tsv
│       ├── md5sum.txt
│       ├── pbta-cnv-consensus-gistic.zip
│       ├── pbta-cnv-consensus.seg.gz
│       ├── pbta-fusion-putative-oncogenic.tsv
│       ├── pbta-histologies.tsv
│       ├── pbta-sv-manta.tsv.gz
│       ├── pbta_germline_plp_calls.tsv
│       ├── release-notes.md
│       ├── snv-consensus-plus-hotspots.maf.tsv.gz
│       └── snv-mutation-tmb-coding.tsv
├── download_data.sh
└── scripts
    ├── install_bioc.r
    └── install_github.r
```

## Code Authors

Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)) and Ryan Corbett ([@rjcorb](https://github.com/rjcorb))

## Contact

For questions, please submit an issue or send an email to Jo Lynne Rokita: rokita@chop.edu

