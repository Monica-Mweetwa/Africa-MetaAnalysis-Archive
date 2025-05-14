# Africa-MetaAnalysis-Archive
This is the archive for code used to generate results reported in the paper titled: A META-ANALYSIS OF GUT MICROBIOME RESEARCH IN MALNOURISHED AFRICAN POPULATIONS: A NATURAL LANGAUGE PROCESSING APPROACH

* The first step was extraction of microbes from json files which was exceuted in python using the code file: Extract_annotations.ipynb
* The next step was summarisation of study charcterisrtics implemented in pythin using the code file: PythonAnnotation_Dec.ipynb
* Phylogenetic tree creation and annotation was done in 2 ways:
* A) Phylogenetic tree of all studies included stratified by region and sequencing method: PhyloTrees_12052025.Rmd
* B) Phylogenetic tree of all studies that compared healthy and undernourished children: CaCo_PhyPlot_12052025.R

The imput data used is stored in the Data folder.
* The data used for phylogenetic tree constraction is stored in the 'Data/All_Studies' folder
* The case control analysis used data stored in the 'Data/CaseControl_Analysis'

