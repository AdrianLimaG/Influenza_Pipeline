# Workflow 3 Script
_______________________________________

## Run Nextclade

<br />

#nextcalde.py

- **nextclade_runner()**
    - Calls combine_nextclade_output, nextclade_data_obj, run_nextclade
    - This function runs nextcalde and pushes results to Influenza DB

<br />

#nextcalde_helper.py

- **run_nextclade()**
    - Runs nextcalde from resources path
    - Returns Dictionary of postive hits

- **concat_fasta()**
    - Concates all fasta files per samples to be passed on to nextcalde

- **combine_nextclade_output()**
    - combines nextcalde output into one file


- **nextclade_data_obj() CLASS**
    - Initialize class
    - Requires path to data directory to import sql db credentials 

    - **get_nextclade_dfs()**
        -Read in combined nextcalde file and create DF

    - **database_push()**
        - Pushes  nextcalde DF to Influenza DB

    - **setup_db()**
        - Sets up Influenza DB connection


<br />