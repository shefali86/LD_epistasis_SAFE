NOTE: SAFE downloaded package to run my analyses. Performed a git clone here to get all information. Please see below to run some analyses


INTRODUCTION
============

SAFE (or Spatial Analysis of Functional Enrichment) is an automated network annotation algorithm. Given a biological network and a set of functional groups or quantitative features of interest, SAFE performs local enrichment analysis to determine which regions of the network are over-represented for each group or feature. SAFE visualizes the network and maps the detected enrichments onto the network.

SAFE was originally implemented in MATLAB and stored at  <https://bitbucket.org/abarysh/safe/>. However, as of early 2017, the MATLAB implementation is only maintained for legacy purposes. New work related to SAFE is moving  to Python and this repository. 

**WARNING. This package is still in development. Please use caution.**

Go on SAFE project to look at the initial guidelines for installation if required. Here are steps I performed:

```
cd safepy-master/
source activate safepy_env_
jupyter-notebook 
import sys
from os.path import expanduser

# Add path to folder containing safepy
sys.path.append('/Users/ssetia/Downloads/safepy-master/')

import safe
%matplotlib inline
import numpy as np

sf = safe.SAFE()

```


My steps on file preparation for SAFE

```

## Biofilter to annotate gene, upstream and downstream
biofilter.py -v --grch-build-version 37 -k ~/group/software/biofilter/2.4.0/loki-b37.db --position-file all_ohta_models_chr_pos.txt -a position gene upstream downstream --overwrite --verbose --prefix all_ohta_models_chr_pos

## Extract closest gene from Biofilter file
cat all_ohta_models_chr_pos2.position.gene-upstream-downstream |awk -F '\t' '{if($4=="" && ($6<$8 || $8==""))print $1":"$3,$5; else if($4=="" && ($6>$8 || $6==""))print $1":"$3,$7; else if($4!="")print $1":"$3,$4}' >all_ohta_models_chr_pos.closestgene.txt


##Map SNP-SNP to Gene-Gene and weight for network edge file where our weight should be combination of BA and CVC

## Add background based on file `/project/ritchie09/personal/shefali_lpc/eMERGE_LD_Ohta_analyses/1KG_data/allpaiwise_models_epistatic_selection_all.txt` where weight is 20 and replace significant weights based on each phenotype in eMERGE or UKB results.

join -1 3 -2 1 <(join -1 1 -2 1 <(cat all_ohta_models_chr_pos.closestgene.txt |sort -k1,1) <(cat all_ohta_models_weight_250000.txt |sort -k1,1) |sort -k3,3) <(cat all_ohta_models_chr_pos.closestgene.txt |sort -k1,1) |awk '{print $3,$5,$4}' >all_ohta_models_chr_pos.closestgene.weigts.txt 
## Note there is a gene name called 3.8-1.5 and it seems like an actual gene and not error


## Since Cytoscape is crashing when I am trying to load this file to make network, I did following to aggregate same gene-gene models

>data<-read.table("all_ohta_models_chr_pos.closestgene.weigts.txt",header=T)
>data2 =aggregate(data, by =list(data$Gene1, data$Gene2), FUN=mean)
> write.table(data2, file="all_ohta_models_closestgene.weights.aggregated.txt",sep=" ", quote=FALSE, row.names=FALSE)
>q()

cat all_ohta_models_closestgene.weights.aggregated.txt |awk '{print $1,$2,$5}' >all_ohta_models_closestgene.weights.aggregated2.txt

## Generate all combination weights. We should not have to do this part. So ignore
join -e '5e-08' -a1 -a2 -o1.1,1.2,2.1 <(cat all_ohta_models_closestgene.weights.aggregated2.txt |awk '{print $1":"$2,$3}' |sort -k1,1) <( sh permutations.sh |awk '{if($1!=$2) print $1":"$2}' |sort -k1,1) |awk '{print $3,$2}' |sed 's/:/ /g' |grep -v Gene >all_ohta_models_closestgene.weights.aggregated2.allcomb.txt

## The above command is running from very long, so instead of waiting I tried one more way: IGNORE
sh permutations.sh >all_ohta_models_closest.allcom.txt

grep -Fwvf <(cat all_ohta_models_closestgene.weights.aggregated2.txt |awk '{print $1,$2}') 
all_ohta_models_closest.allcom.txt |awk '{print $1,$2,"5e-08"}' >all_ohta_models_closest.allcomb.remaining.lowestvalue.txt

cat all_ohta_models_closestgene.weights.aggregated2.txt all_ohta_models_closest.allcomb.remaining.lowestvalue.txt |grep -v Gene > all_ohta_models_closestgene.weights.aggregated2.allcomb.try2.txt

## START HERE AGAIN Generate Node attribute file (in yeast, they map genes to ORF, I am just leaving this file as same gene-gene names)
cat all_ohta_models_chr_pos.closestgene.txt |awk '{print $2}' |sort -u |awk '{print $1,$1}' >all_ohta_models_chr_pos.closestgene.node.attribute.txt

## Map Genes to GO Terms and then generate attribute file
biofilter.py -v --grch-build-version 37 -k ~/group/software/biofilter/2.4.0/loki-b37.db --gene-file all_ohta_models.closest.genes.txt -a gene group source --source go --overwrite --verbose --prefix all_ohta_models.closest.genes.GO

cat all_ohta_models.closest.genes.GO.gene.group-source|sed 's/ /_/g' |awk -F '\t' '{if($3=="go") print $1,$2,"1"}' |sort -u >all_ohta_models.closest.genes.GOterms.txt


>data<-read.table ("all_ohta_models.closest.genes.GOterms.txt",header=T)
> library(reshape)
> data2<-cast(data, ORF~GO, fill=0)
Using VALUE as value column.  Use the value argument to cast to override this choice
> write.table(data2, file="all_ohta_models.closest.genes.GOterms_matrix.txt",sep="\t", quote=FALSE, row.names=FALSE)
>q()

## Steps to re-create files where we only include Gnees that map to GO terms
cat all_ohta_models.closest.genes.GOterms_matrix.txt |awk '{print $1}' |sort -u >Genes_withGO.txt

grep -Fwf Genes_withGO.txt all_ohta_models_chr_pos.closestgene.txt >all_ohta_models_chr_pos.closestgene_wGO.txt

join -1 3 -2 1 <(join -1 1 -2 1 <(cat all_ohta_models_chr_pos.closestgene_wGO.txt |sort -k1,1) <(cat all_ohta_models_weight_250000.txt |sort -k1,1) |sort -k3,3) <(cat all_ohta_models_chr_pos.closestgene_wGO.txt |sort -k1,1) |awk '{print $3,$5,$4}' >all_ohta_models_chr_pos.closestgene.weights_wGO.txt

>data<-read.table("all_ohta_models_chr_pos.closestgene.weights_wGO.txt",header=T)
>data2 =aggregate(data, by =list(data$Gene1, data$Gene2), FUN=mean)
>write.table(data2, file="all_ohta_models_closestgene.weights_wGO.aggregated.txt",sep=" ", quote=FALSE, row.names=FALSE)
>q()

cat all_ohta_models_closestgene.weights_wGO.aggregated.txt |awk '{print $1,$2,$5}' >all_ohta_models_closestgene.weights_wGO.aggregated2.txt

sh permutations.sh >all_ohta_models_closest.wGO.allcom.txt

grep -Fwvf <(cat all_ohta_models_closestgene.weights_wGO.aggregated2.txt |awk '{print $1,$2}') all_ohta_models_closest.wGO.allcom.txt |awk '{print $1,$2,"5e-08"}' >all_ohta_models_chr_pos.closestgene.wgo.remaining.lowestvalue.txt

cat all_ohta_models_closestgene.weights_wGO.aggregated2.txt all_ohta_models_chr_pos.closestgene.wgo.remaining.lowestvalue.txt |grep -v Gene |sed 's/ /\t/g' >all_ohta_models.closestgene.wgo.allpairwise.aggregated.txt

```

HELP
====

Please direct all questions/comments regarding SAFE to Anastasia Baryshnikova (<abaryshnikova@calicolabs.com>).

The main repository for this code is at <https://github.com/baryshnikova-lab/safepy>. Please subscribe to the repository to receive live updates about new code releases and bug reports.


For this git page direct questions/comments to Shefali S. Verma (shefali.setiaverma@pennmedicine.upenn.edu).


CITE SAFE
==========

The manuscript describing SAFE and its applications is available at:

> Baryshnikova, A. (2016). Systematic Functional Annotation and Visualization of Biological Networks. Cell Systems. <http://doi.org/10.1016/j.cels.2016.04.014>
