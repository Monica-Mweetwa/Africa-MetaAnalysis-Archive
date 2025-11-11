# Africa-MetaAnalysis-Archive
This is the archive for code used to generate results reported in the paper titled: A META-ANALYSIS OF GUT MICROBIOME RESEARCH IN MALNOURISHED AFRICAN POPULATIONS: A NATURAL LANGAUGE PROCESSING APPROACH

HTML files from pubmed central were donwloaded and converted to json files using [AutoCorpus](https://github.com/omicsNLP/Auto-CORPus). The microbe names were extracted using [MicrobELP](https://github.com/omicsNLP/microbELP).
* The first step at annotation was extraction of microbes from json files which was executed in python using the code file: Extract_annotations.ipynb
* The next step was summarisation of study characteristics manually extracted from text which was implemented in python using the code file: PythonAnnotation_Dec.ipynb
* Phylogenetic tree creation and annotation was done in 3 ways:
* A) Phylogenetic tree of all studies included in the analysis: PhyloTrees_12052025.Rmd
* B) Phylogenetic tree of all studies included stratified by region and sequencing method: PhyloTrees_12052025.Rmd
* C) Phylogenetic tree of all studies that compared healthy and undernourished children: CaCo_PhyPlot_12052025.R

The imput data used is stored in the 'Data' folder.
* The data used for phylogenetic tree constraction is stored in the 'Data/All_Studies' folder
* The case control analysis used data stored in the 'Data/CaseControl_Analysis'

