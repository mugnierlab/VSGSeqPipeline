To convert these fasta files to blast databases:  

```
makeblastdb -in vsgs_427_cds_030414.fa -input_type fasta -dbtype nucl -out tb427_vsgs
```
```
makeblastdb -in NOTVSGs_wflanks_070617.fa -input_type fasta -dbtype nucl -out NOTvsgs
```
```
makeblastdb -in 3570cds_with_duplicates_annotated.fa -input_type fasta -dbtype nucl -out EATRO1125_vsgs
```
```
cat vsgs_427_cds_030414.fa 3570cds_with_duplicates_annotated.fa > concatEATRO1125_tb427VSGs.fa 
```
```
makeblastdb -in concatEATRO1125_tb427VSGs.fa -input_type fasta -dbtype nucl -out concatAnTattb427
```
