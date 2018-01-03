04/21/2017 - Downloaded GENIE PNET clinical data from: http://www.cbioportal.org/genie
  - Cancer type: Pancreatic cancer #463
  - Cancer Type detailed: Pancreatic neuroendocrine #76
    - File name: ~/Desktop/PughLab/NET-seq/datasets/genie/genie_data/GENIE_pnet_clinicalData.txt
    - 76 samples / 72 patients
  - MEN1, DAXX, ATRX patients:
    - File name: ~/Desktop/PughLab/NET-seq/datasets/genie/genie_data/GENIE_pnet_daxx-atrx-men1_clinicalData.txt
    - 39 samples / 38 patients
    
04/21/2017 - Downloaded GENIE CNV data from SYNAPSE: https://www.synapse.org/#!Synapse:syn7851253
  - File name: ~/Desktop/PughLab/NET-seq/datasets/genie/genie_data/data_CNA.txt
 
04/21/2017 - Downloaded GENIE mutation data from SYNAPSE: https://www.synapse.org/#!Synapse:syn7851252
  - File name: ~/Desktop/PughLab/NET-seq/datasets/genie/genie_data/data_mutations_extended.txt
  
04/21/2017 - Downloaded GENIE CNV seg file data from SYNAPSE: https://www.synapse.org/#!Synapse:syn7851253.23
  - File name: ~/Desktop/PughLab/NET-seq/datasets/genie/genie_data/genie_public_segments.seg

04/21/2017 - Used grep to cut the sample names from each clinicalData sheet to pull out seg and mutation information.

head -1 genie_public_segments.seg > GENIE_pnet_cna.seg
rm GENIE_pnet_mutations.txt
for i in $(cut -f5 GENIE_pnet_clinicalData.txt); do
  grep $i genie_public_segments.seg >> GENIE_pnet_cna.seg
  grep $i data_mutations_extended.txt >> GENIE_pnet_mutations.txt
done

#  grep $i genie_public_segments.seg >> GENIE_pnet_cna_genes.txt

head -1 genie_public_segments.seg > GENIE_pnet_DAM_cna.seg
rm GENIE_pnet_DAM_mutations.txt
for i in $(cut -f5 GENIE_pnet_daxx-atrx-men1_clinicalData.txt); do
  grep $i genie_public_segments.seg >> GENIE_pnet_DAM_cna.seg
  grep $i data_mutations_extended.txt >> GENIE_pnet_DAM_mutations.txt
done