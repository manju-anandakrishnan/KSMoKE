# This repository maintains the code base for the web interface described in the manuscript, 'Improving kinase enrichment analysis with predicted associations and an interactive kinase centric network explorer'.

The application can be accessed at https://research.bioinformatics.udel.edu/KSMoFinder/ <br>

**1. Conduct kinase enrichment analysis and explore connections between enriched kinase(s) and substrate proteins** <br/><br/>
<img src=images/readme_figure2.png width="800">

<br/> 

**2. Retrieve predicted kinases of a phosphorylation site and explore biological connections between the kinase and substrate protein** <br/>
<img src=images/readme_figure1.png width="800">

# Application set up
**Create a conda environment with the required libraries** <br/>
```
conda create env --name netvisKE --file requirements.txt
```

**Activate the conda environment** <br/>
```
conda activate netvisKE
```

**Update .env file with db configurations** <br/>

**Launch the application** <br/>
```
streamlit run home.py
```





