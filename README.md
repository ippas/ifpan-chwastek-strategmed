# ifpan-chwastek-strategmed
RNA-seq of osteoarthritis rat model

1. Samples were downloaded from Novogene file transfer system
2. file list was generated : `ls data/*/*/*/*fq.gz > fq-files.txt`
3. QC was run in batches of 25 samples:
```
cat fq-files.txt | head -25 | xargs -I {} -n 1 docker run --rm -d -v $PWD:/data pegi3s/fastqc /data/{}
```
4. To run the multiQC report:
```
docker run -d --rm -v $PWD:/data ewels/multiqc:latest multiqc /data -o /data
```
5. The rest of the analysis was run with the (IntelliSeq workflow)[https://gitlab.com/intelliseq/workflows/raw/rna-seq-paired-end@1.10.0/src/main/wdl/pipelines/rna-seq-paired-end/latest/rna-seq-paired-end.wdl]

Get the input:
```
files_1=`ls data/*/*/*/*_1.fq.gz | xargs -i echo \"{}\",`
files_2=`ls data/*/*/*/*_2.fq.gz | xargs -i echo \"{}\",`
samples=`ls data/*/*/*/*_1.fq.gz | xargs -i bash -c 'BASENAME=$(echo {} | cut -d "." -f 1 | cut -d "/" -f 4); echo \"$BASENAME\",'`

echo "{\"rna_seq_paired_end.fastqs_left\":[$files_1],\"rna_seq_paired_end.fastqs_right\":[$files_2],\"rna_seq_paired_end.sample_names\":[$samples]}" > input.json
```
don't forget the extra commas and add this text:
```
    "rna_seq_paired_end.organism_name": "Rattus norvegicus",
    "rna_seq_paired_end.release_version": "101",
    "rna_seq_paired_end.analysis_id": "rna_seq_paired_end_test"
}
```

5. Analysis was run in R 3.4 (see analysis code here)[]
