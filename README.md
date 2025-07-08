# Germline pathogenic variation impacts somatic alterations and patient outcomes in pediatric CNS tumors

Ryan J. Corbett^, Rebecca Kaufman^, Shelly W. McQuaid, Zalman Vaksman, Saksham Phul, Miguel A. Brown, Jennifer L. Mason, Sebastian M. Waszak, Bo Zhang, Chuwei Zhong, Emily Blauel, Heena Desai, Ryan Hausler, Ammar S. Naqvi, Jessica M. Daggett, Alex Sickler, Evan C. Cresswell-Clay, Patricia J. Sullivan, Antonia Chroni, Zhuangzhuang Geng, Elizabeth M. Gonzalez, Yuankun Zhu, Allison P. Heath, Marilyn Li, Penn Medicine BioBank, Regeneron Genetics Center, Phillip B. Storm, Adam C. Resnick, Kara N. Maxwell, Kristina A. Cole, Angela J. Waanders, Miriam Bornhorst, Suzanne MacFarland, Jo Lynne Rokita+, Sharon J. Diskin+

^Equal authorship

+Co-senior authors

:tada: This work is now preprinted on [MedRxiv](https://www.medrxiv.org/content/10.1101/2025.02.04.25321499v1). :newspaper:


__DISCLAIMER__: This repository contains de-identified data only, and does NOT contain any raw germline variant data.  

## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct analyses for the manuscript noted above.

1. Clone the repository
```
git clone git@github.com:diskin-lab-chop/pbta-germline-somatic.git
```

2. Pull the docker container:
```
docker pull pgc-images.sbgenomics.com/diskin-lab/pbta-germline:1.0.1
```

3. Start the docker container

From the `pbta-germline-somatic` folder, run:

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/pbta-germline-somatic pgc-images.sbgenomics.com/diskin-lab/pbta-germline:1.0.0
```
 
Users can also run Rstudio in the project docker container from a web browser using the instructions below:

__Local Development in Rstudio__ (Max OS X and Linux users only)

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/pbta-germline-somatic pgc-images.sbgenomics.com/diskin-lab/pbta-germline:1.0.1
```

Then, navigate to `localhost:8787` in your web browser. The username for login is `rstudio` and the password will be whatever password is set in the `docker run` command above (default: `pass`).

__Development using Amazon EC2, depending on your open ports__

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 80:8787 -v $PWD:/home/rstudio/pbta-germline-somatic pgc-images.sbgenomics.com/diskin-lab/pbta-germline:1.0.1
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
│   ├── v11
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

7. To run all analyses provided in the manuscript (needs to be run on an HPC or EC2):
```
bash scripts/run_analysis_modules.sh
```

## Code Authors

Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)) and Ryan Corbett ([@rjcorb](https://github.com/rjcorb))

## Contact

For questions, please submit an issue.

