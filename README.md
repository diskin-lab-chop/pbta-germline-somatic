# Germline variation across OpenPBTA patients

## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct analyses for the manuscript noted above.

1. Clone the repository
```
git clone git@github.com:diskin-lab-chop/OpenPBTA-germline.git
```

2. Pull the docker container:
```
docker pull jrokita1/pbta-germline:version0.1
```

3. Start the docker container, from the `OpenPBTA-germline` folder, run:
```
docker run --name container_name -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/OpenPBTA-germline jrokita1/pbta-germline:version0.1
```

4. To execute shell within the docker image, from the `OpenPBTA-germline` folder, run:
```
docker exec -ti container_name bash
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
│   └── collapse-tumor-histologies
│       ├── 01-collapse-tumor-histologies.R
│       └── results
│           ├── primary-tumor-histologies-collapsed-by-germline.tsv
│           └── tumor-histologies-collapsed-by-germline.tsv
├── data
│   └── v1
│       ├── consensus_seg_annotated_cn_autosomes.tsv.gz
│       ├── consensus_seg_annotated_cn_x_and_y.tsv.gz
│       ├── consensus_seg_with_status.tsv
│       ├── independent-specimens.rnaseq.primary-plus-polya.tsv
│       ├── independent-specimens.rnaseq.primary-plus-stranded.tsv
│       ├── independent-specimens.wgs.primary-plus.tsv
│       ├── independent-specimens.wgs.primary.tsv
│       ├── independent-specimens.wgswxs.primary-plus.tsv
│       ├── independent-specimens.wgswxs.primary.tsv
│       ├── md5sum.txt
│       ├── pbta-cnv-consensus-gistic.zip
│       ├── pbta-cnv-consensus.seg.gz
│       ├── pbta-fusion-putative-oncogenic.tsv
│       ├── pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
│       ├── pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
│       ├── pbta-germline-all-20220729.tsv
│       ├── pbta-histologies.tsv
│       ├── pbta-sv-manta.tsv.gz
│       ├── release-notes.md
│       ├── snv-consensus-plus-hotspots.maf.tsv.gz
│       └── snv-mutation-tmb-coding.tsv
├── download_data.sh
├── figures
└── scripts
    ├── install_bioc.r
    └── install_github.r
```

## Code Authors

Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)) and Ryan Corbett ([@rjcorb]((https://github.com/rjcorb))

## Contact

For questions, please submit an issue or send an email to Jo Lynne Rokita: rokita@chop.edu

