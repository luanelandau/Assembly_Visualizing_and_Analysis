# Assembly_Visualizing_and_Analysis
For these scripts, I am mainly visualizing the draft assembly as dot-plots against other assemblies. It is useful to have an idea of how big are your contigs, the coverage you have, or if you are interested in seeing how a specific part of your assembly looks like.

DotPlots.sh : If you wish to create dotplots against a reference genome for the entire assembly or maybe just one chromosome. 
Example:


<img width="613" alt="image" src="https://github.com/user-attachments/assets/26acdcd8-bfa9-4abf-8a5e-ce404a50c4c7">


dotplots_hap1_hap2_targetregion.sh : this is for optimizing dotplots for regions of interest. It was specifically created to look at the amylase gene copies (SVs) in human genomes. But can be used for any region of interest. 
Example:


<img width="540" alt="image" src="https://github.com/user-attachments/assets/cde295af-7884-4186-91f5-513666219aae">


align_region_of_interest.sh : dotplots for a region of interest. It is the previous version of the above, but it could be useful if you dont need to optimize the dotplots (e.g. only one contig aligns perfectly with your region). You might want to try this one before trying the above.
extract_aligned_fasta_from_assembly.sh : If you want to extract a fasta file from your assembly matching a specific region of a reference genome. Say you want to look at a region that contains genes or repeats. 

