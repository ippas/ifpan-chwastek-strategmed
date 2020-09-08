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
5. Align with Hisat2 to the rat genome:

Get te inputs:
```
ls data/*/*/*/*_1.fq.gz | xargs -i bash -c 'BASENAME=$(echo {} | cut -d "." -f 1 | cut -d "_" -f 1,2,3); echo $BASENAME' | xargs -i bash -c 'echo "{\"align_to_rat_genome.align_with_hisat2.fastq1\":\"{}_1.fq.gz\",\"align_to_rat_genome.align_with_hisat2.sample_name\":\"{}\",\"align_to_rat_genome.align_with_hisat2.fastq2\":\"{}_2.fq.gz\"}">{}-input.json'
```

Run two samples at the time: 
```
ls data/*/*/*/*_1.fq.gz | head -2 | xargs -i bash -c 'BASENAME=$(echo {} | cut -d "." -f 1 | cut -d "_" -f 1,2,3); echo $BASENAME' | xargs -i bash -c 'java -jar /opt/tools/cromwell-44.jar run https://raw.githubusercontent.com/gosborcz/workflows/master/align-with-hisat2-to-rat-genome -i {}-input.json > log-{}.txt'
```
