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
