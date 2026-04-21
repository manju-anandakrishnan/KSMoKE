# This repository maintains the code base for the web interface described in the manuscript, KSMoKE.

The application can be accessed at https://research.bioinformatics.udel.edu/KSMoFinder/ <br>

**1. Conduct kinase enrichment analysis and explore connections between enriched kinase(s) and substrate proteins** <br/><br/>
<img src=images/readme_figure2.png width="800">

<br/> 

**2. Retrieve predicted kinases of substrate using KSMoFinder and explore biological connections between kinase and substrate protein** <br/>
<img src=images/readme_figure1.png width="800">

# Application set up
**Create a conda environment with the required libraries** <br/>
```
conda create env --name ksmoke --file requirements.txt <br/>
```


**Activate the conda environment** <br/>
```
conda activate ksmoke <br/>
```

**Update .env file with db configurations** <br/>

**Launch the application** <br/>
```
streamlit run home.py
```





