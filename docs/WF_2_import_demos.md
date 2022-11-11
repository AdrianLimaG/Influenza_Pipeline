# Workflow 2 Script
_______________________________________

## Import Patient Demographic Information

<br />

#import_demo.py

- **run_import_demo()**
    - Calls demographical import class and runs
    - Require Dictionary to paths, and path to resources folder, and iram output dir
    - Creates json file to store assembly metrics, coverage and depth


<br />

#import_demo_helper.py

- **demographics_import() CLASS**
    - Initialize class
    - Requires path to data directory to import both sql db credentials and samples_metrics
    - Creates a DF from samples metrics json file 

- **get_lims_demographics()**
    - Imports patient demographical information from HORIZON
    - Import PCR results from HORIZON
    - Save demographical information to results folder for laster use
    - 

- **format_lims_df()**
    - Formats HORIZON for Influenza DB

- **merge_dfs()**
    - Umerges both lims_df with sample metrics df

- **database_push()**
    - Pushes merged DF to Influenza DB

- **setup_db()**
    - Sets up Influenza DB connection


<br />
