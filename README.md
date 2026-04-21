# This repository maintains the code base for the web interface described in the manuscript, KSMoKE.

The application can be accessed at https://research.bioinformatics.udel.edu/KSMoFinder/ <br>

**Conduct kinase enrichment analysis and explore connections between enriched kinase(s) and substrate proteins** <br/>
<img src=images/readme_figure2.png width="800">

**Create a conda environment with the required libraries** <br/>
conda create env --name ksmoke --file requirements.txt <br/>


**Activate the conda environment** <br/>
conda activate ksmoke <br/>

**Update .env file with db configurations** <br/>

**Launch the application** <br/>
streamlit run home.py






