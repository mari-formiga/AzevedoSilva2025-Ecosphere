Description of the data and file structure
Between November 2012 and April 2013, woody plants at seven sites of Cerrado sensu stricto were searched for ant-plant interactions. In each locality, 5 transects of 200 m were established, at least 1 km apart from one another. Transects were divided into 20 plots of 10 m length and data were collected in 10 alternate plots. Each transect was characterized by its climatic, vegetational and resource variables. For Single Nucleotide Polymorphisms (SNPs) genotyping, we selected *Camponotus crassus *specimens from samples collected during fieldwork. For data collection details, please see Material and Methods.

Files and variables
File: Ecological_data-reviewed.txt
Description: Local, ID, and coordinates of transects (used as sample unit for analyses). 
For each transect is also available the operational explanatory variables. 
Variables:
Transect = transect ID
Local = sampling site. Bra = Brasília; Can = Canastra; Cip = Cipó; Iti = Itirapina; Ser = Serra Azul; Ema = Emas; Vea = Veadeiros
Lat = transect latitude (unit: decimal degree; WGS94)
Long = transect longitude (unit: decimal degree; WGS94)
D.mean = Mean vegetation density (unit: mean number of plants/m2)
prop.res= Proportion of plants with sugar-rich resources (extrafloral nectaries and/or trophobionts)
T= Mean annual temperature (unit: C x 10)
P = Annual precipitation (unit: mm)
S.ants = ants rerified species richness
prop.crassus = proportion of plants with Camponotus crassus

File: Genetic_data (VCF file)
Description: 186 SNPs for 143 individuals of C. crassus.

File: code-reviewed.R
Description: Code used to run the analyses described in Material and Methods.


