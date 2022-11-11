# Workflow 0 Script
_______________________________________

## Formats Raw Fastq Files from MinION for down stream analysis

<br />

#merge_fastq.py

- **merge_seq_fastq()**
    - Calls find_samples from helper file
    - Path to Fastq files and path to samplesheet required

<br />

#merger_fastq_helper.py

- **find_samples()**
    - Function organises barcode folders produced from MinION
    - Matches barcodes with SampleID from sample sheet
    - Calls read_sample_sheet and concat_fastq to help
    -Returns Dictionary of all{SampleID:Path_to_combined_fastq}

- **concat_fastq()**
    - Concats the multiple fastq files produced per sample into one fastq file per samples
    - Returns {SampleID:Path_to_combined_fastq} for given samples

- **read_sample_sheet()**
    - Reads in sample sheet
    - Returns dictionary of {barcode:sampleID}



<br />
