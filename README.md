# SpliceAlnViz
Version 2019

Interactive alignment visualization tool
highlight important components of gene  structures such as exons, introns, splice sites, and predicted exons.​
----------

Authors:  Safa Jammali
Université de Sherbrooke, Canada
Cobius Lab:  https://cobius.usherbrooke.ca/
for questions email us at Safa.Jammali@USherbrooke.ca


### Requirements:


-Bio
-skBio
-pandas
-xlwt
-Ete toolkit
-Blast+ standalone
-argparse


###Dash Installation (https://dash.plot.ly/installation?_ga=2.50536041.1535164263.1551149387-2095754676.1539110927)

pip install dash==0.38.0  # The core dash backend
pip install dash-html-components==0.13.5  # HTML components
pip install dash-core-components==0.43.1  # Supercharged components
pip install dash-table==3.5.0  # Interactive DataTable component 


### Running SpliceAlnViz:
python2.7 MSAViz.py

###Usage

MSA:
1. Upload input files
2. Click OK
3. Choose features
4. Choose sequences
5. Write height and width value

Search motif:
1. Write motif sequence
2. Click submit
3. download the results

Artificial MSA:
1. Write height and width value
2. Select All/CDS/Gene

 
