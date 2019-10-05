# Scripts used in Van Paemel et al. (2019)

## Introduction
These are the scripts that were used in the manuscript "Low-cost, non-invasive classification of pediatric solid tumors using reduced representation bisulfite sequencing of cell-free DNA: a proof-of-principle study". See also materials and method section of the manuscript for more information.

The classifier was built using publicly available data. Unfortunately, due to data-transfer agreements, it's not possible to cannot share all the source files necessary to reproduce all the results and you will have to request access to the original authors. 

## Publicly available datasets
The datasets can be obtained from:
- Healthy cfDNA: [EGAS00001000566](https://www.ebi.ac.uk/ega/studies/EGAS00001000566)  & [GSE79277](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79277)
- White blood cell: [E-MTAB-4187](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-4187/)
- Rhabdomyosarcoma: [EGAS00001000884](https://www.ebi.ac.uk/ega/studies/EGAS00001000884)
- Ewing sarcoma: [GSE89041](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89041), [GSE97529](https://www.ncbi.nlm.nih.gov/gds/?term=GSE97529[Accession]), [GSE90496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90496)
- Malignant rhabdoid tumor: [GSE107946](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107946)
- Intracranial tumors: [GSE90496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90496)
- Neuroblastoma: TARGET-L3: ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/NBL/methylation_array/
- Osteosarcoma: TARGET-L3: ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/OS/methylation_array/
- Wilms tumor: TARGET-L3: ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/WT/methylation_array/
- Clear cell sarcoma of the kidney: TARGET-L3: ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/CCSK/methylation_array/

## Preprocessing raw sequencing reads
After the sequencing run has finished, all fastq.gz files are put into one folder and the snakemake pipeline in the preprocessing folder (`code/preprocessing/Bismark_pipeline_SE.snakefile`) was initiated.
- For RRBS data: reads were not deduplicated according to the Bismark user guide.
- For WGBS data: reads were deduplicated.

The .cov.gz files from running bismark methylation extractor were then used either to classify according to the histopathological diagnosis (cfRRBS data) or to add to the reference set (WGBS data).

## Get MspI regions from GRCh37 with mkrrgenome
A .bed file with MspI regions between 20 and 200 bp in the human genome is needed for the steps below. Such a .bed file can be made with [mkrrgenome](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378906/).

This is possible with the `-m` flag (sizes are specified with `-M`, default 40-220):

```bash
./mkrrgenome -g ./GRCh37-lite_chr/ -M 20,200 -H ./GRCh37-lite_chr/rr_ -m 20,200 > RRBS_regions.txt
```

However, the output of this is not in `.bed` format, but can be made with `mkrrgenome2bed.py`.

```bash
python mkrrgenome2bed.py -i [input] -o [output]
```

## Get the intersection between HM450K and RRBS with MakeNewClusters.ipynb
CpGs were clustered in order to use more mappable reads. We adjusted the target regions to make them more suited for RRBS data. We merged all remaining regions within 1 bp from each other with BEDtools. Finally, clusters were retained if they contain at least 3 CpGs covered on the Illumina HM450K array, resulting in 14103 clusters covering 61750 probes on the HM450K array (`RRBS_450k_intersectClusters.tsv`, included in this repository).

The script requires the Infinium HM450K manifest file, HumanMethylation450_15017482_v1-2.csv, and be download from [the Illumina website](http://emea.support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html).

## Make reference set with MakeTrainTest.py
`MakeTrainTest.py` builds the input files for `runMeth_atlas.py` (both the reference set and the samples that need to be classified). It takes the argument `--train` for rebuilding the reference set and `--viz` for tSNE and UMAP plots (`--viz` only works if `--train` is enabled). Normally, the reference set should only be built once (but every time a sample or tumor entity is added, it needs to be updated). 

The output files (in the example train and test_methy and test_depth) have the same structure as described in the [CancerLocator repository](https://github.com/jasminezhoulab/CancerLocator).

Note:
- The reference set is referred to as "train" in the scripts, but no 'training' is happening like in traditional machine learning algorithms.
- Because the publicly available data all comes in different formats (sometimes a m x n matrix, sometimes m x 1 matrix, sometimes something completely different), every tumor entity has it's own way of importing. There is an example file for each tumor entity to illustrate what the original input file looked like. The beta values in these example files were scrambled.

## NNLS wrapper: deconvolve_resids.py
The wrapper around the NNLS algorithm was cloned from the meth atlas paper repository (https://github.com/nloyfer/meth_atlas), no adjustments to this script were made.

## Classification with runMeth_atlas.py 
After running the previous steps, there should be at least two files present:

- trainFile: Each line represents a reference sample. The first column is the sample type and the remaining columns are methylation values (beta values) of the clusters.
- testBetaFile: Each line represents a test sample. The first column is the sample type and the remaining columns are methylation values (beta values) of the clusters.

```python
# Classification for plasma samples
inputFile_tumor_atlas = "./resources/train"
inputFile_samples = "./resources/test_beta"
```

The reference sets used in the manuscripts can be found under `code/train_methatlas_*_manuscript.csv`.

After running the script, the output can be found under `./classifySamples/output/classification/classificationResults_methAtlas.csv`

## Copy number aberrations
WisecondorX (https://github.com/CenterForMedicalGeneticsGhent/WisecondorX) was used to generate copy number aberration profiles for both sWGS and cfRRBS data. For this, mapping vs. GRCh38 (as opposed to GRCh37 in the other sections) was done. The reason for not mapping everything vs. GRCh38 was mainly due to compatibility issues, because the Illumina HM450K platform is completely GRCh37-based.
### cfRRBS
Steps taken to obtain cfRRBS CNA profiles:
1. Convert .bam files (after mapping with Bismark and sorting with samtools) to .npz files. with 200kb bins: `WisecondorX convert --binsize 200000 ${BAM} ${OUTPUT_DIR}/${SAMPLE}.npz`
2. A reference set with healthy controls samples was made. It's essential that this reference set is done with the same DNA isolation, library prep and sequencing machine to account for systematic errors: `WisecondorX newref reference_input_dir/*.npz reference_output.npz --nipt --binsize 400000`
3. Make CNA calls: 
```
for npz in *_bt2.npz; 
  do 
    SAMPLE=${npz%%.*}
    WisecondorX predict $npz reference.GRCh38.400kb.npz ./output400kb/$SAMPLE --beta 0.15 --minrefbins 50 --plot --bed
  done
```
### sWGS
The steps required to obtain sWGS CNA profiles and IchorCNA tumor fraction estimations are described in more detail in `code/preprocessing/sWGS_pipeline_SE.snakefile`. 

## Downstream analysis
See RMarkdown notebook for more details. 
