# BTI summer intern

## Name: 

## Date: June-27-2023

# Table of Contents
- [BTI summer intern](#bti-summer-intern)
  - [Name:](#name)
  - [Date: June-27-2023](#date-june-27-2023)
- [Table of Contents](#table-of-contents)
  - [Install miniconda and others](#install-miniconda-and-others)
  - [Trim adaptor and low-quality sequences, polyA and rRNA to collect only the good reads](#trim-adaptor-and-low-quality-sequences-polya-and-rrna-to-collect-only-the-good-reads)
  - [Mapping good reads to reference genome](#mapping-good-reads-to-reference-genome)
    - [Map reads to Pennellii](#map-reads-to-pennellii)
    - [Map reads to M82](#map-reads-to-m82)
  - [Compare mismatch counts for each read aligned to the M82 and Pennellii genome to assign each read to the genome with a lower mismatch value](#compare-mismatch-counts-for-each-read-aligned-to-the-m82-and-pennellii-genome-to-assign-each-read-to-the-genome-with-a-lower-mismatch-value)
    - [Python script to average all lines with the same IDs](#python-script-to-average-all-lines-with-the-same-ids)
      - [Standard script](#standard-script)
      - [Using multi-threading](#using-multi-threading)
      - [Multi processing](#multi-processing)
    - [Assign reads to the genome for which they have a lower average mismatch count](#assign-reads-to-the-genome-for-which-they-have-a-lower-average-mismatch-count)
  - [Find number of ids mapped to Pennellii and M82](#find-number-of-ids-mapped-to-pennellii-and-m82)
  - [Extract reads mapped to Pennellii and M82 based on ids mapped to each parent using mismatch counts](#extract-reads-mapped-to-pennellii-and-m82-based-on-ids-mapped-to-each-parent-using-mismatch-counts)
    - [Uses gatk package to extract corresponding reads](#uses-gatk-package-to-extract-corresponding-reads)
    - [Python script to extract corresponding reads](#python-script-to-extract-corresponding-reads)
      - [Script without multi-processing](#script-without-multi-processing)
      - [Python script using multi processing](#python-script-using-multi-processing)
    - [Extracts and maps reads using python program](#extracts-and-maps-reads-using-python-program)
    - [Counts number of reads mapped to M82 and Pennellii](#counts-number-of-reads-mapped-to-m82-and-pennellii)
  - [Assign mapped reads from each parent to genes using FeatureCounts](#assign-mapped-reads-from-each-parent-to-genes-using-featurecounts)
  - [Collect metrics for normalization and extract fragment length from metrics](#collect-metrics-for-normalization-and-extract-fragment-length-from-metrics)
    - [Collect metrics](#collect-metrics)
    - [Extract fragment length](#extract-fragment-length)
  - [Liftoff](#liftoff)
    - [Code to remove genes with multiple pairings](#code-to-remove-genes-with-multiple-pairings)



## Install miniconda and others
The following code to be run from the command line installs the necessary packages for cleaning the reads from the F1 of the M82 and Pennellii tomatoes. 

```bash
# Working directory
cd /data/mukund/tools

# Download miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --no-check-certificate

# Install miniconda
bash ./Miniconda3-latest-Linux-x86_64.sh

# Export to environment
export PATH="/data/mukund/tools/miniconda3/bin:$PATH"

# install hisat2
conda install -c conda-forge mamba
mamba install -c bioconda hisat2

# install trimmomatic
mamba install -c bioconda trimmomatic

# install bowtie
mamba install -c bioconda bowtie

```

##  Trim adaptor and low-quality sequences, polyA and rRNA to collect only the good reads 
This script removes adaptor, low-quality, polyA, and rRNA from the reads. The trimmomatic package is used to clean the adaptor and low-quality sequences, PRINSEQ++ is used to remove the polyA sequences, and bowtie is used to remove the rRNA sequences, since only mRNA sequences are needed for analysis. It also includes code to return the number of clean reads following the cleaning process. 

```bash
# Working directionary
cd /data/mukund/M82_pennellii/01trim

# Get all the sample ids
ll *R1.fastq.gz | awk '{print $9}' | sed 's/_R1.fastq.gz//g' > sample.id
# Process the data
for i in $(cat ../00fastq/sample.id); do
    # Trim adaptor and low quality
    trimmomatic PE ../00fastq/${i}_R1.fastq.gz ../00fastq/${i}_R2.fastq.gz \
        ${i}_R1P.fq ${i}_R1U.fq \
        ${i}_R2P.fq ${i}_R2U.fq \
        -threads 60 \
        ILLUMINACLIP:/data/mukund/tools/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
        SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40
    # Remove polyA
    prinseq++ -threads 60 -VERBOSE 0 \
        -fastq ${i}_R1P.fq \
        -fastq2 ${i}_R2P.fq \
        -min_len 40 \
        -trim_tail_left 10 \
        -trim_tail_right 10 \
        -out_name ${i}
    # Remove rRNA
    bowtie -v 3 -k 1 -p 60 --al ${i}_R1P.rRNA.fq \
        /data/zhaojiantao/database/rRNA_combined/rRNA_combined \
        -q ${i}_good_out_R1.fastq >/dev/null
    bowtie -v 3 -k 1 -p 60 --al ${i}_R2P.rRNA.fq \
        /data/zhaojiantao/database/rRNA_combined/rRNA_combined \
        -q ${i}_good_out_R2.fastq >/dev/null
    grep @ ${i}_R1P.rRNA.fq | cut -f1 -d ' ' | sed 's/@//g' > ${i}_R1P.rRNA.readsID 
    grep @ ${i}_R2P.rRNA.fq | cut -f1 -d ' ' | sed 's/@//g' > ${i}_R2P.rRNA.readsID
    cat ${i}_R1P.rRNA.readsID ${i}_R2P.rRNA.readsID | sort | uniq > ${i}_rRNA.sorted.readsID
    wc -l ${i}_rRNA.sorted.readsID
    seqkit grep -vf ${i}_rRNA.sorted.readsID ${i}_good_out_R1.fastq -o ${i}_R1P.clean.fq.gz
    seqkit grep -vf ${i}_rRNA.sorted.readsID ${i}_good_out_R2.fastq -o ${i}_R2P.clean.fq.gz
    # Remove temporary files
    rm *tq *ID *fq 
done

# Save the above command line to Qualityclean.sh and run: nohup sh Qualityclean.sh > Qualityclean.log &
# Extract the summary Statistics
# Total raw reads
grep "Input Read Pairs" Qualityclean.log | awk '{print $4}'
# After trimming
grep "Input Read Pairs" Qualityclean.log | awk '{print $7}'
# Count the number of PolyA
for i in $(cat sample.id); do
        echo $(cat ${i}_single_out_R1.fastq ${i}_bad_out_R1.fastq|wc -l)/4|bc
done

# Cross-check with the number of remained reads

for i in $(cat sample.id); do
        echo $(cat ${i}_good_out_R1.fastq|wc -l)/4|bc
done

# Count the number of cleaned reads
for i in $(cat sample.id); do
        echo $(zcat ${i}_R1P.clean.fq.gz|wc -l)/4|bc
done

# Remove unnecessary files
rm *R1P.fq *R2P.fq *rRNA.fq *readsID *single* *bad* *good* *R1U.fq *R2U.fq

```

## Mapping good reads to reference genome
The following code provides a pipeline to map the clean reads to the reference genome of the M82 and Pennellii tomatoes. It will be necessary to change the file paths of the reference genomes to the proper file paths. It also provides code to calculate the number of reads mapped to each genome. 

Here, sample.id is a list of the Pennellii sample names. 
### Map reads to Pennellii

```bash
# Working directory
cd /data/zhaojiantao/Carmen/Carmen_ASE/Pennellii/02map2Pennellii

# Index the reference genome
hisat2-build -f /data/zhaojiantao/Carmen/References/Pennellii/Pennellii.fa \
        /data/zhaojiantao/Carmen/References/Pennellii/Pennellii.fa

# Map the reads
for i in $(cat ../01trim/sample.id); do
        #Map the reads
        hisat2 --rna-strandness RF -p 30 --no-softclip \
                --summary-file $i.clean.summary \
                -x /data/zhaojiantao/Carmen/References/Pennellii/Pennellii.fa \
                -1 ../01trim/${i}_R1P.clean.fq.gz \
                -2 ../01trim/${i}_R2P.clean.fq.gz \
                -S $i.clean.sam
        # Convert sam to bam and extract no mismatches
        samtools view -@ 100 -S -b $i.clean.sam > $i.clean.bam
        # Remove sam file
        rm $i.clean.sam
done

# Total clean reads
for i in $(cat Pennellii.list); do
        sed -n '1p' $i.clean.summary | awk '{print $1}'
done

# Calculate the total reads mapped to the reference genome
for i in $(cat Pennellii.list); do
        awk '{print $1}' $i.clean.summary | tr -s "\n" "\t" | awk -v OFMT='%.0f' '{print $4+$5+$8+$13/2+$14/2}'
done

```


### Map reads to M82
The two code snippets marked "same function as..." perform the same function and return the total clean reads. It is not necessary to run both. Note that hisat and samtools are both running on 100 threads in this code snippet, which may need to be changed depending on the capabilities of your system. 

Here, M82.list is a list of the M82 samples. 
```bash
# Working directory
cd /data/zhaojiantao/Carmen/Carmen_ASE/Pennellii/02map2M82

# Index the reference genome
hisat2-build -f /data/zhaojiantao/Carmen/References/M82/M82.fa \
        /data/zhaojiantao/Carmen/References/M82/M82.fa

# Map the reads
for i in $(cat M82.list); do
        #Map the reads
        hisat2 --rna-strandness RF -p 100 --no-softclip \
                --summary-file $i.clean.summary \
                -x /data/zhaojiantao/Carmen/References/M82/M82.fa \
                -1 ../01trim/${i}_R1P.clean.fq.gz \
                -2 ../01trim/${i}_R2P.clean.fq.gz \
                -S $i.clean.sam
        # Convert sam to bam and extract no mismatches
        samtools view -@ 100 -S -b $i.clean.sam | samtools sort -@ 100 > $i.clean.sort.bam
        # Remove sam file
        rm $i.clean.sam
done

# Total clean reads - same function as 2 below
for i in $(cat M82.list); do
        sed -n '1p' $i.clean.summary | awk '{print $1}'
done

# Calculate the total reads mapped to the reference genome
for i in $(cat M82.list); do
        awk '{print $1}' $i.clean.summary | tr -s "\n" "\t" | awk -v OFMT='%.0f' '{print $4+$5+$8+$13/2+$14/2}'
done

# Count total reads - same function as 2 above
for i in $(cat M82.list.complete); do
        cut -d " " -f 1 $i.clean.summary | head -1 
done
```


## Compare mismatch counts for each read aligned to the M82 and Pennellii genome to assign each read to the genome with a lower mismatch value
For each read, extracts and compares the number of base pairs that are mismatched with each genome and assigns the read to the genome with lower mismatch count.

### Python script to average all lines with the same IDs
Averages the mismatch counts for the same IDs in each genome for later comparison and assignment. Multi-threading and multi-processing capabilities are provided, although multi-threading may not work as intended due to the possibility of chunks being cut between the same ids. The multi-processing script will prevent this, and is highly recommended for efficiency.

#### Standard script
Uses a single thread. 
```py
import sys

# Trigger to prevent running of program if proper arguments are not provided
run = True
# Collect arguments from command line for input file and output file
arguments = sys.argv

# Provide help if called with unexpected carguments
if (arguments[1] == "h" or arguments[1] == "help"):
        print("This script takes arguments of the form python mismatch.compare.py <input> <output>. The input file must be a sorted, tab-delimited list of read IDs and their corresponding mismatch counts. It will average the mismatch counts for each read and write each read, along with its average mismatch count, to the output file.")
        run = False
elif (len(arguments) != 3):
        print("This script takes arguments of the form python mismatch.compare.py <input> <output>. The input file must be a sorted, tab-delimited list of read IDs and their corresponding mismatch counts. It will average the mismatch counts for each read and write each read, along with its average mismatch count, to the output file.")
        run = False
if (run):
        # Open provided files 
        with open(arguments[1], "r") as data, open(arguments[2], "w") as output:
                # Read lines from data file and collect first id
                lines = data.readlines()
                lastid = lines[0].split('\t')[0]
                total = 0
                count = 0
                # Iterate through lines of input file
                for line in lines:
                        # Collect current id
                        line = line.strip()
                        id = line.split('\t')[0]
                        id = id.strip()

                        # Collect current mismatch
                        mismatch = line.split('\t')[1].strip()

                        # Continue adding mismatch counts until id does not match previous, then write average and corresponding lastid to output file
                        if (id == lastid):
                                count = count + 1
                                total = total + int(mismatch)
                        else:
                                output.write(lastid + "\t" + str(total/count) + "\n")
                                count = 1
                                total = int(mismatch)
                        lastid = id
```

#### Using multi-threading
Multi-threading is not recommended and may provide incorrect output. For capability to run on multiple CPUs, use multi-processing script instead.
```py
import sys
import threading

# Trigger to prevent running of program if proper arguments are not provided
run = True

# Collect arguments from command line for input file and output file
arguments = sys.argv

# Provide help if called with unexpected arguments
if (arguments[1] == "h" or arguments[1] == "help"):
    print("This script takes arguments of the form python mismatch.compare.py <input> <output>. The input file must be a sorted, tab-delimited list of read IDs and their corresponding mismatch counts. It will average the mismatch counts for each read and write each read, along with its average mismatch count, to the output file.")
    run = False
elif (len(arguments) != 3):
    print("This script takes arguments of the form python mismatch.compare.py <input> <output>. The input file must be a sorted, tab-delimited list of read IDs and their corresponding mismatch counts. It will average the mismatch counts for each read and write each read, along with its average mismatch count, to the output file.")
    run = False

if run:
    # Define averager function to average lines with the same id
    def averager(lines):
        # Collect first id and initialize counter and sum variables
        lastid = lines[0].split('\t')[0]
        total = 0
        count = 0

        # Iterate through provided lines
        for line in lines:
            # Collect current line and id
            line = line.strip()
            id = line.split('\t')[0]
            id = id.strip()
            # Collect current mismatch
            mismatch = line.split('\t')[1].strip()
            # Compare current id with previous id
            # If equal, then continue incrementing total and count variables
            # Else write last id and corresponding average to output file
            if id == lastid:
                count = count + 1
                total = total + int(mismatch)
            else:
                output.write(lastid + "\t" + str(total / count) + "\n")
                count = 1
                total = int(mismatch)
            lastid = id
    
    # Open provided files
    with open(arguments[1], "r") as data, open(arguments[2], "w") as output:
        # Read lines from provided input file
        lines = data.readlines()
        # Set number of threads to use
        num_threads = 4  
        # Determine size of chunks to be averaged on each thread
        chunk_size = len(lines) // num_threads
        threads = []
        # Determine chunks
        # May cut chunks between two equal ids and not provide complete average
        # Use multi-processing for accurate result
        for i in range(num_threads):
            start = i * chunk_size
            end = (i + 1) * chunk_size if i < num_threads - 1 else None
            chunk = lines[start:end]
            t = threading.Thread(target=averager, args=(chunk,))
            threads.append(t)
            t.start()

        # Join threads
        for t in threads:
            t.join()
```

#### Multi processing
Will provide the most efficient method of averaging mismatch values. Changing the variable num_processes will change the number of required CPUs. 
```py
import sys
import multiprocessing

# Trigger to prevent running of program if proper arguments are not provided
run = True

# Collect arguments for input file and output file from command line
arguments = sys.argv

# Provide help if called with unexpected arguments
if (arguments[1] == "h" or arguments[1] == "help"):
    print("This script takes arguments of the form python mismatch.compare.py <input> <output>. The input file must be a sorted, tab-delimited list of read IDs and their corresponding mismatch counts. It will average the mismatch counts for each read and write each read, along with its average mismatch count, to the output file.")
    run = False
elif (len(arguments) != 3):
    print("This script takes arguments of the form python mismatch.compare.py <input> <output>. The input file must be a sorted, tab-delimited list of read IDs and their corresponding mismatch counts. It will average the mismatch counts for each read and write each read, along with its average mismatch count, to the output file.")
    run = False

if run:
    # Define function to be called in each process to average lines in a chunk
    def averager(chunk):
        # Collect first id and intialize total and count variables for computing averages
        lastid = chunk[0].split('\t')[0]
        total = 0
        count = 0
        output_lines = []
        # Iterate over lines in each chunk
        for line in chunk:
            # Collect current line and id
            line = line.strip()
            id = line.split('\t')[0]
            id = id.strip()
            # Collect current mismatch
            mismatch = line.split('\t')[1].strip()
            # If current id is equal to the last id, then continue incrementing total by the mismatch count and increasing the counter variable by 1
            # Else write last id and corresponding average to the output file 
            if id == lastid:
                count += 1
                total += int(mismatch)
            else:
                output_lines.append(lastid + "\t" + str(total / count))
                count = 1
                total = int(mismatch)
            lastid = id
        # Append very last id to ensure that it is included
        output_lines.append(lastid + "\t" + str(total / count))
        return output_lines

    # Open provided data
    with open(arguments[1], "r") as data:
        # Read lines from provided data 
        lines = data.readlines()
        # Number of processes to use 
        num_processes = 10
        # Compute size of chunks
        chunk_size = len(lines) // num_processes
        chunks = []
        start = 0
        # Create chunks
        # Only ends chunks when first id of new chunk does not match lastid of previous chunk to ensure all ids are averaged
        for i in range(num_processes):
            end = None
            if i < num_processes - 1:
                for j in range(start + chunk_size, len(lines)):
                    if lines[j].split('\t')[0] != lines[j-1].split('\t')[0]:
                        end = j
                        break
            chunk = lines[start:end]
            chunks.append(chunk)
            start = end or start + chunk_size
    # Create multi processing pool and run averager function with each chunk
    pool = multiprocessing.Pool(processes=num_processes)
    results = pool.map(averager, chunks)
    pool.close()
    pool.join()

    # Write output from all chunks to output file 
    with open(arguments[2], "w") as output:
        for result_chunk in results:
            output.write('\n'.join(result_chunk) + '\n')
```

### Assign reads to the genome for which they have a lower average mismatch count
Extracts only the ids and corresponding mismatch counts from each file. Then, for each id, uses the provided multi-processing script above to average the mismatch counts for each parent, compare the mismatch counts for each parent, and assigns the read to the parent with a lower mismatch count. Note that the averaging script is called multi.mismatch.compare.py here, but the program name should be changed to the name of your corresponding program. Also note that the sametools command is using 20 threads, but this may need to be changed depending on the capabilities and limitations of your processing system. 

Here, sample.id is a list of all sample names. 
```bash
for i in $(cat sample.id); do

        # Read bam files for Pennellii mapping and create a file of ids and corresponding mismatch counts
        samtools view -@ 20 ../02map2Pennellii/$i.clean.bam | grep "NM:i" | cut -f 1 > $i.clean.bam.ids.Pennellii
        samtools view -@ 20 ../02map2Pennellii/$i.clean.bam | grep -o "NM:i.." | sed "s/NM:i://g" > $i.clean.bam.mismatch.indicators.Pennellii
        paste $i.clean.bam.ids.Pennellii $i.clean.bam.mismatch.indicators.Pennellii | sort --parallel=30 > $i.clean.bam.mismatches.Pennellii

        # Read bam files for M82 mapping and create a file of ids and corresponding mismatch counts
        samtools view -@ 20 ../02map2M82/$i.clean.sort.bam | grep "NM:i" | cut -f 1 > $i.clean.bam.ids.M82
        samtools view -@ 20  ../02map2M82/$i.clean.sort.bam | grep -o "NM:i.." | sed "s/NM:i://g" > $i.clean.bam.mismatch.indicators.M82
        paste $i.clean.bam.ids.M82 $i.clean.bam.mismatch.indicators.M82 | sort --parallel=30 > $i.clean.bam.mismatches.M82

        # Average mismatch counts for reads with the same id using multi-processing script 
        python multi.mismatch.compare.py $i.clean.bam.mismatches.Pennellii $i.clean.bam.mismatches.Pennellii.average
        python multi.mismatch.compare.py $i.clean.bam.mismatches.M82  $i.clean.bam.mismatches.M82.average

        # Join two files into one file with id number, average for Pennelli on left, and average for M82 on the right
        join $i.clean.bam.mismatches.Pennellii.average $i.clean.bam.mismatches.M82.average > $i.clean.bam.combined

        # Compare mappings for each ID and create two files with ids assigned to each parent
        awk -v output1="$i.ids.map2Pennellii" -v output2="$i.ids.map2M82" '{ if ($2 < $3) print $1 > output1; else if ($3 < $2) print $1 > output2 }' "$i.clean.bam.combined"

        # Remove temporary files
        rm $i.clean.bam.ids.Pennellii $i.clean.bam.mismatch.indicators.Pennellii $i.clean.bam.mismatches.Pennellii $i.clean.bam.ids.M82 $i.clean.bam.mismatch.indicators.M82 $i.clean.bam.mismatches.M82 $i.clean.bam.mismatches.Pennellii.average $i.clean.bam.mismatches.M82.average $i.clean.bam.combined
        done
```

## Find number of ids mapped to Pennellii and M82 
Extracts the number of read ids mapped to each parent. Here, sample.id is a list of all sample names.
```bash
for i in $(cat sample.id); do 
        wc -l $i.ids.map2M82 | cut -d ' ' -f 1

for i in $(cat sample.id); do
        wc -l $i.ids.map2Pennellii | cut -d ' ' -f 1 

```
## Extract reads mapped to Pennellii and M82 based on ids mapped to each parent using mismatch counts
Extracts reads from original bam file and assigns them to each parents based on the assignment of their corresponding id. Two methods are provided, one using the package gatk, and the other using a python script. The python script is slightly more efficient and produces smaller files, but the content of the output is the same. The python script may require a virtual environment, and does require the BAM files to be sorted and indexed, however, which may increase computing time if not already done. 

### Uses gatk package to extract corresponding reads
Runs on one thread. sample4.id is the list of samples. 
```bash
# Creates two files with the bam information of mapped ids 
for i in $(cat sample4.id); do
        /data/zhaojiantao/tools/gatk-4.2.6.1/gatk FilterSamReads \
        -I ../02map2M82/$i.clean.sort.bam \
        -O $i.map2M82.bam \
        --READ_LIST_FILE $i.ids.map2M82 \
        --FILTER includeReadList
done
```

### Python script to extract corresponding reads

#### Script without multi-processing
Does not run more efficiently than gatk, multiple threads are necessary to achieve superior computing time. 
```py
# Requires virtual environment: conda activate myenv

import pysam
import sys

# Read arguments from command line
arguments = sys.argv

def extract_reads_by_ids(bam_file, ids_file, output_bam):
    # Open BAM files
    in_bam = pysam.AlignmentFile(bam_file, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", header=in_bam.header)

    # Read ids from file
    with open(ids_file, "r") as ids:
        id_set = set(line.strip() for line in ids)

    # Extract reads matching the ids
    for read in in_bam.fetch():
        if read.query_name in id_set:
            out_bam.write(read)

    # Close BAM files
    in_bam.close()
    out_bam.close()

# Call function using command line arguments
extract_reads_by_ids(arguments[1], arguments[2], arguments[3])
```


#### Python script using multi processing
Note that a higher number of processes will not necessarily reduce the computing time and may even increase it. Depending on the size of the file, a specific number of processes greater than one ususally provides peak efficiency. 
```py
import pysam
import sys
import multiprocessing

# Define function to extract ids
def extract_reads_by_ids(bam_file, ids_file, output_bam):
    # Open BAM files
    in_bam = pysam.AlignmentFile(bam_file, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", header=in_bam.header)

    # Read ids from file
    with open(ids_file, "r") as ids:
        id_set = set(line.strip() for line in ids)

    # Extract reads matching the ids
    for read in in_bam.fetch():
        if read.query_name in id_set:
            out_bam.write(read)

    # Close BAM files
    in_bam.close()
    out_bam.close()

if __name__ == '__main__':
    # Collect arguments from command line
    arguments = sys.argv
    # Create multi-processing pool
    pool = multiprocessing.Pool()
    # Define number of processes to be used
    numprocesses = 30
    # Call function with multi-processing
    pool.starmap(extract_reads_by_ids, [(arguments[1], arguments[2], arguments[3]) for i in range(numprocesses)])
    pool.close()
    pool.join()

```


### Extracts and maps reads using python program
The following code uses the provided multi-processing python script to extract the reads mapped to Pennellii and M82. Note that they will need to be sorted and indexed first, and the number of threads used for this can be changed in the samtools sort and samtools index parameters. 
```bash
# Prepare bam files for python code by sorting and indexing
for i in $(cat sample.id); do
        samtools sort -@ 20 $i.clean.bam > $i.clean.sorted.bam
        samtools index -@ 20 $i.clean.sorted.bam
        rm $i.clean.bam
done

# Run python script to extract reads mapped to Pennellii and M82 and write them to new bam files
for i in $(cat sample.id); do
        python multiextraction.py ../02map2M82/$i.clean.sort.bam $i.ids.map2M82 $i.map2M82.bam
        python multiextraction.py ../02map2Pennellii/$i.clean.sorted.bam $i.ids.map2Pennellii $i.map2Pennellii.bam
done

``` 

### Counts number of reads mapped to M82 and Pennellii
```bash
#M82
for i in $(cat sample.id); do 
        samtools flagstat -@ 30 $i.map2M82.bam | head -7 | tail -n 1 | cut -d ' ' -f 1
        done
# Pennellii
for i in $(cat sample.id); do 
        samtools flagstat -@ 30 $i.map2Pennellii.bam | head -7 | tail -n 1 | cut -d ' ' -f 1
        done


```

## Assign mapped reads from each parent to genes using FeatureCounts
```bash
# Use FeatureCounts to map M82 reads to genes
/data/zhaojiantao/tools/subread-2.0.3-Linux-x86_64/bin/featureCounts -T 64 -s 2 -p --countReadPairs -M \
        -a /data/zhaojiantao/Carmen/Carmen_ASE/Pennellii/04orthologs/M82.uniq_specific.gff3 \
        -o RawReads2M82.count \
        -t gene \
        -g ID \
        -O \
        ../03map2Parents/*map2M82.bam \
        2> RawReads2M82.count.log

# Use FeatureCounts to map Pennellii reads to genes
/data/zhaojiantao/tools/subread-2.0.3-Linux-x86_64/bin/featureCounts -T 64 -s 2 -p --countReadPairs -M \
        -a /data/zhaojiantao/Carmen/Carmen_ASE/Pennellii/04orthologs//Pennellii.uniq_specific.gff3 \
        -o RawReads2Pennellii.count \
        -t gene \
        -g ID \
        -O \
        ../03map2Parents/*map2Pennellii.bam  \
        2> RawReads2Pennellii.count.log 
```

## Collect metrics for normalization and extract fragment length from metrics

### Collect metrics
```bash
# Collect metrics from M82 
for i in $(cat sample.id) ; do 
java -jar /data/zhaojiantao/tools/picard/build/libs/picard.jar  CollectInsertSizeMetrics  I=$i.map2M82.bam         O=$i.map2M82.metrics.txt         M=0.5         H=test.pdf & done

# Collect metrics from Pennellii
for i in $(cat sample.id) ; do 
java -jar /data/zhaojiantao/tools/picard/build/libs/picard.jar  CollectInsertSizeMetrics  I=$i.map2Pennellii.bam         O=$i.map2Pennellii.metrics.txt         M=0.5         H=test.pdf & done
```

### Extract fragment length 
```bash
# Pennellii
for i in $(cat sample.id); do 
        head -8 $i.map2Pennellii.metrics.txt | tail -1 | cut -f 1
done

# M82 
for i in $(cat sample.id); do 
        head -8 $i.map2M82.metrics.txt | tail -1 | cut -f 1
done
```



## Liftoff
Note that this does not relate to the read analysis above, and can be used to remove genes with multiple pairings after using the Liftoff program. 
### Code to remove genes with multiple pairings
```py
with open("PennelliiAllLiftedData", "r") as alldata, open("Pennelliicleanedliftedata", "w") as cleandata:
    alllines = alldata.readlines()
    backid = ""
    nowid = ""
    nextid = ""
    lastline = ""
    # Compare each id with the one before and after it, and only include it if it does not match either
    for line in alllines:
        line = line.strip()
        nextid = line.split('\t')[0]
        nextid = nextid.strip()
        if (nowid != nextid and nowid != backid):
            cleandata.write(lastline + "\n")
        backid = nowid
        nowid = nextid
        lastline = line
    
    cleandata.write(lastline)
```

[def]: #bti-summer-intern
