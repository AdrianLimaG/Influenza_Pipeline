
# General Package for State Lab Influenza A Clade Detection
_______________________________________

## The package contains the following workflows in their respective subdirectories:

<br />

### **Workflow 0:** [Merge Sequencing Data](docs/WF_0_merge_seq_data.md)
 - Joins all individual samples fastq files into one fastq file per sample
 - Map barcodes to sample IDs
 - Create/Return sample dictionary holding path to fastq 

  > This step is required, since the assembler only takes in 1 fastq file per sample.<br>
  > And need to map samples to their respective sample ID.

<br />

<br />

### **Workflow 1:** [Assemble Reads using IRMA](docs/WF_1_irma.md)
 - Run [IRMA assembler](https://wonder.cdc.gov/amd/flu/irma/)
 - After consense sequences are created, check protein assembly stats and write them to a JSON file 

  > This step is required, to obatin Influenza consense sequences.<br>

<br />


<br />

### **Workflow 2:** [Import demographics](docs/WF_2_import_demos.md)
 - Open HORIZON LIMS database (Oracle).
 - Join all demographics with sample ID.
 - Join demographics with assembly stats.
 - Push new demographics to Influenza DB MS SQL database.
 - Write demographical iformation for final result file

  > This step is absolutely required, since the sample ID is the primary key in the database<br>
  > making it impossible to insert any other results further down the workflow.

<br />

<br />

### **Workflow 3:** [Run Nextclade](docs/WF_3_nextclade.md)
 - Run Nextcalde to determine clade of seqeunced influenza virus
 - Check against Influenza A h3n2,h1n1pdm types (can be modified to include FLU B)
 - Parse the nextcalde data
 - Push the Nextclade data to Influenza DB MS SQL database.
   
<br />
<br />

### **Workflow 4:** [Gisaid Report](docs/WF_4_gisaid_export.md)
 - Currently NOT IMPLETMENTED DUE TO GISAID FLU being down.
 - Extract required data from nextclade files.
 
  > Current not implemented due to GISAID Flu DB being down
  
<br />
<br />

### **Workflow 5:** [Build Epi Report](docs/WF_5_final_report.md)
 - Ask user for search to perform
 - Takes both Nextclade results and demographic information and builds report 
 - Formats consense fasta files
 - Aligns samples to referance Influenza A h3n2,h1n1pdm type based on nextclade restuls
 - Builds phylogenetic tree for each influenza A types based on alignement files
 - Cleans ups intermediate files


<br />

