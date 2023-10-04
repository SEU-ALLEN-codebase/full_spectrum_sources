# Local Morphology Generation (Optional)

The first step in constructing the microenvironment involves the generation of morphologies, which can encompass dendrites, complete morphologies, or any sub-neuronal arbors derived from full morphologies. In this project, we have utilized automatically traced local morphologies to represent local image blocks. As such, the initial phase necessitates the acquisition of an adequate number of morphologies from neuronal images. Below, we provide a comprehensive set of exemplary procedures for tracing neurons within extensive neuronal blocks.

## Automatic Tracing Codes Location

The automatic tracing codes can be found at `generation/42k_script`. Prior to utilizing these codes, please ensure that the following prerequisites are met:

- Cropped images with dimensions of 512x512x256 voxels (in xyz order).

## Procedure

1. **Removal of Dense Images with Multiple Somas**
   
   - Location: `generate/42k_script/rm_multi`
   - This is an optional step and is only necessary if there are numerous image blocks with multiple somas.
   - To perform this step, you should provide at least one text file containing the coordinates of all somas. Examples of such files include "/home/vkzohj/data/enhanced_soma_images/all.txt" or the pickled file "coords.pickle".
   - You can estimate the grid-based soma counts matrix using the script "statis.py".
   - Finally, using these files, you can filter out all dense image blocks and save the remaining image blocks into a text file "passed.txt" using a similar approach as demonstrated in the script "filter.py".

2. **Enhance the Image Using "imPreProcess" Plugin**
   
   - Location: `generation/42k_script/enh/run.sh`
   - The enhancement method "imPreProcess" is bound with Vaa3D.
   - This enhancement process is detailed in the paper "Image enhancement to leverage the 3D morphological reconstruction of single-cell neurons" (Guo et al., Bioinformatics, 2022, [doi: 10.1093/bioinformatics/btab638](https://doi.org/10.1093/bioinformatics/btab638)).

3. **Tracing with APP2**
   
   - Location: `generation/42k_script/app2`
   - APP2 is a popular automatic tracing algorithm (Xiao et al., bioinformatics, 2013, [doi: 10.1093/bioinformatics/btt170](https://doi.org/10.1093/bioinformatics/btt170)), and it is available in Vaa3D.
   - Exemplary batch running codes can be found in the specified directory.

4. **Tracing with neuTube**
   
   - Location: `generation/42k_script/ntb`
   - Similar to step 3 (APP2 tracing), but using another popular auto-tracing algorithm neuTube (Feng et al., eNeuro, [doi: 10.1523/ENEURO.0049-14.2014](https://doi.org/10.1523/ENEURO.0049-14.2014)).
   - Batch running codes are located in the specified directory.

5. **Pruning the APP2 Reconstructions**
   
   - Location: `generation/42k_script/seg_prune`
   - Refer to the script "run.sh" for guidance on pruning APP2 reconstructions.

6. **Consensus Calculation**
   
   - Location: `generation/42k_script/consensus`
   - Estimation of the consensus morphology between pruned APP2 reconstruction and neuTube reconstruction.
   - Relevant codes can be found in the specified directory.




# Microenvironment construction




