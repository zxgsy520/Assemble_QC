# Assemble_QC
Eukaryotic Genome Assembly Evaluation Software

## Statistics of assembled genomes
```
./stat_genome genome.fasta -o stat_genome.tsv
./stat_genome genome.fasta -o stat_genome.tsv --split >contig.fasta
```
#StatType	ContigLength	ContigNumber	ScaffoldLength	ScaffoldNumber	GapLength	GapNumber
N50	12,925	257	44,702	89	400	521
N60	9,893	378	33,996	125	200	641
N70	7,404	538	27,567	170	100	784
N80	5,172	757	19,033	231	400	917
N90	3,277	1,084	11,444	321	200	1,042
Longest	178,915	1	281,293	1	500	1,121
Total	13,667,795	1,789	13,780,305	663	112,510	1,126
Length>=1kb	13,664,395	1,779	13,780,305	663	112,510	1,126
Length>=2kb	13,105,799	1,394	13,602,323	528	112,410	1,125
Length>=5kb	11,040,997	777	13,295,289	434	111,610	1,117
![image](https://user-images.githubusercontent.com/36355222/190681636-5e82e801-3043-4597-b46f-263059aa3823.png)

