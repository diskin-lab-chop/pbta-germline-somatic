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

3. Start the docker container, from the `pbta-germline-somatic-integration` folder, run:
```
docker run --name container_name -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/OpenPBTA-germline jrokita1/pbta-germline:version0.1
```

4. To execute shell within the docker image, from the `PBTA-ALT-analysis` folder, run:
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
│   └── generate-histologies
├── data
├── download_data.sh
├── figures
└── scripts
    ├── install_bioc.r
    └── install_github.r
```

## Code Authors

Jo Lynne Rokita ([@jharenza](https://github.com/jharenza))

## Contact

For questions, please submit an issue or send an email to Jo Lynne Rokita: rokita@chop.edu

