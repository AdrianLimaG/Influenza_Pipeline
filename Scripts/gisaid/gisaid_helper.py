from ms_sql_handler import ms_sql_handler
from ..reader import read_txt
import os
import datetime
import pandas as pd
import reader






class gisaid_obj():
    # constructor
    def __init__(self,cache_path):
        demo_cahce= reader.read_json(cache_path+"/data/gisaid.json")
        for item in [*demo_cahce] :
            setattr(self,item, demo_cahce[item])

    def setup_db(self):
        self.db_handler = ms_sql_handler(self)
        self.db_handler.establish_db()


    def get_priority(self):
        if not self.reportable:
            print("\n==========================================",
                        "\nALERT! You are generating a report for",
                        "\n------------------"
                        "\nSURVEILLANCE ONLY."
                        "\n------------------",
                        "\nIf you wish to report the results to",
                        "\nHORIZON, please change the 'reportable' variable",
                        "\nin the static_cache.json file to a 1.",
                        "\n==========================================\n")
        print("\nGetting list of priority samples...")
        self.setup_db()
        samples = self.db_handler.ss_read(query=self.read_query_tbl1_priority).values.astype(str).tolist()
        self.priority_lst = [sample[0].strip() for sample in samples]
        print(" Done!\n")


    def scan_db(self):
        # in this function, we need to get the maximum gisaid number
        # so we know how to label future samples
        # we also will scan the database for the last seven days for
        # samples eligible for upload based on the qc cutoffs set
        # in data/static_cache.json
        self.setup_db()
        self.next_gisaid = int(self.db_handler.ss_read(query=self.read_query_tbl1_max_gisaid).iat[0, 0]) + 1
        prev_week = (datetime.date.today() - datetime.timedelta(days = 31)).strftime("%Y%m%d")
        self.read_query_tbl1_eligible_hsn = self.read_query_tbl1_eligible_hsn.replace("{prev_week}", prev_week)
        if not self.reportable:
            self.read_query_tbl1_eligible_hsn = self.read_query_tbl1_eligible_hsn.replace("{reportable}", "in (0, 1)")
        else:
            self.read_query_tbl1_eligible_hsn = self.read_query_tbl1_eligible_hsn.replace("{reportable}", " = " + str(self.reportable))
        self.hsn_lst = self.db_handler.sub_read(query=self.read_query_tbl1_eligible_hsn)['hsn'].astype(str).tolist()
        if len(self.hsn_lst) == 0:
            raise ValueError("==================================================================================\nError:\nNo eligible samples for gisaid report!! - All samples already reported\n==================================================================================\n")

    def get_gisaid_df(self):
        # use the hsn list from above to get all hsns into single dataframe
        self.gisaid_start = self.db_handler.sub_lst_read(query=self.read_query_tbl1, lst=self.hsn_lst)
        # self.user now stored in the private_cache.json file
        #self.user = input("\nPlease input the user for this report\n--> ")
        self.hsn_dict = {'hsn':[], 'gisaid':[]}

    def compile_fasta(self, run_id):
        # create the name of the file that will hold the completed fasta data
        # make both destination files for metadata and fasta
        self.file_no = 1
        date1 = datetime.datetime.strptime(run_id[7:17], "%Y-%m-%d").strftime("%m%d%y")
        self.folderpath = self.folderpathbase + "/" + date1 + "/"
        if not os.path.exists(self.folderpath):
            os.makedirs(self.folderpath)
        

        if os.path.exists(self.folderpathbase+"/"+ date1+".txt"):
            print("=========Creating new file number========")
            gisaid_mem = open(self.folderpathbase+"/"+ date1+".txt", "r+")
            #now read the file and find the largers number and set it self.file_no
            self.file_no=len(gisaid_mem.readlines())+1
            #then write this number to mem file
            gisaid_mem.write(str(self.file_no)+"\n")
            gisaid_mem.close()
        else:
            gisaid_mem = open(self.folderpathbase+"/"+ date1+".txt", "w+")
            gisaid_mem.write(str(self.file_no)+"\n")
            gisaid_mem.close()

        self.filepath = date1 + "_" + str(self.file_no) + ".fasta"    
        # get list of fasta files
        self.file_lst = self.gisaid_start["path_to_fasta"].values.tolist()


    def compile_gisaid(self):
        # compile the gisaid template
        self.gisaid_df = pd.DataFrame(self.gisaid_start.rename(columns=self.rename_gisaid_cols_lst))
        self.gisaid_df.sort_values(['wgs_run_date', 'hsn'], inplace=True)
        # self.gisaid_df = add_cols(obj = self,
        #     df = self.gisaid_df,
        #     col_lst = self.add_col_lst,
        #     col_func_map = self.col_func_map)
        self.gisaid_df["hsn"] = self.gisaid_df.apply(lambda row: row['hsn'], axis=1)
        self.gisaid_df.insert(0, "submitter", self.user)
        self.gisaid_df.insert(0, "fn", self.filepath)
        self.gisaid_df["covv_virus_name"] = self.gisaid_df.apply(lambda row: self.get_virus_name(row), axis=1)
        self.gisaid_df.insert(0, "covv_type", "betacoronavirus")
        self.gisaid_df.insert(0, "covv_passage", "Original")
        self.gisaid_df["covv_collection_date"] = self.gisaid_df.apply(lambda row: get_collection_date(row), axis=1)
        self.gisaid_df["covv_location"] = self.gisaid_df.apply(lambda row: self.get_location(row), axis=1)
        self.gisaid_df.insert(0, "covv_add_location", "unknown")
        self.gisaid_df.insert(0, "covv_host", "Human")
        self.gisaid_df.insert(0, "covv_add_host_info", "unknown")
        #This uses self.priorities
        #self.gisaid_df['covv_sampling_strategy'] = self.gisaid_df.apply(lambda row: self.get_comment(row), axis=1)
        self.gisaid_df["covv_gender"] = self.gisaid_df.apply(lambda row: get_sex(row), axis=1)
        self.gisaid_df.insert(0, "covv_patient_status", "unknown")
        self.gisaid_df.insert(0, "covv_specimen", "unknown")
        self.gisaid_df.insert(0, "covv_outbreak", "unknown")
        self.gisaid_df.insert(0, "covv_last_vaccinated", "unknown")
        self.gisaid_df.insert(0, "covv_treatment", "unknown")
        self.gisaid_df.insert(0, "covv_seq_technology", "nanopore minIon")
        self.gisaid_df.insert(0, "covv_assembly_method", self.seq_tech) 
        self.gisaid_df.insert(0, "covv_orig_lab", self.lab_name)
        self.gisaid_df.insert(0, "covv_orig_lab_addr", self.lab_addr)
        self.gisaid_df.insert(0, "covv_provider_sample_id", "unknown")
        self.gisaid_df.insert(0, "covv_subm_lab", self.lab_name)
        self.gisaid_df.insert(0, "covv_subm_lab_addr", self.lab_addr)
        self.gisaid_df.insert(0, "covv_subm_sample_id", "unknown")
        self.gisaid_df.insert(0, "covv_authors", self.authors)
        self.gisaid_df.insert(0, "comment_type", None)
        self.gisaid_df.insert(0, "covv_comment", None)
        # order columns/remove unnecessary columns
        self.gisaid_df = self.gisaid_df[self.full_gisaid_cols_lst]

    def make_fasta_file(self):
        # make the fasta file
        print("\nBuilding the all.fasta file...\n")
        s = ""
        f = open(self.folderpath + self.filepath, "w+")
        for fasta in self.file_lst:
            curr_file = open(fasta, "r")
            file_contents = curr_file.readlines()
            curr_file.close()
            s += "\n\n"
            for line in file_contents:
                if line[0] != ">":
                    s += line
                else:
                    seq_hsn = line[1:8]
                    index = self.gisaid_df[self.gisaid_df['hsn'] == int(seq_hsn)].index.values.tolist()
                    index = index[0]
                    v_name = ">" + self.gisaid_df.iloc[int(index)]['covv_virus_name'] + "\n"
                    s += v_name
            f.write(s)
            s = ""
        f.close()

    def make_gisaid_file(self, run_id):
        date2 = datetime.datetime.strptime(run_id[7:17], "%Y-%m-%d").strftime("%m%d%y")
        templatefilepath = date2 + "_" + str(self.file_no) + "_sql.xlsx"
        self.gisaid_df.to_excel(self.folderpath + templatefilepath, index=False, header=True)

    def database_push(self):
        # put the updated gisaid_nums with the correct samples
        gisaid_df_update = pd.DataFrame.from_dict(self.hsn_dict)
        gisaid_df_update.insert(2, "priority_spec", "0")
        #gisaid_df_update['priority_spec'] = gisaid_df_update.apply(lambda row: get_priority(row, self.priority_lst), axis=1)
        gisaid_df_update_lst = gisaid_df_update.values.astype(str).tolist()
        self.db_handler.lst_ptr_push(df_lst=gisaid_df_update_lst, query=self.write_query_tbl1)

    def get_virus_name(self, row):
        year = str(datetime.datetime.strptime(str(row["doc"]), "%Y-%m-%d").year)
        gisaid = int(self.next_gisaid)
        gisaid_str = f'{gisaid:04}'
        self.hsn_dict['hsn'].append(row["hsn"])
        self.hsn_dict['gisaid'].append(gisaid)
        self.next_gisaid += 1
        return "hCoV-19/USA/" + self.state_abbrev[str(row['state'])] + "-KHEL-" + gisaid_str + "/" + year if str(row['state']).lower() != "unknown" else "hCoV-19/USA/KS" + "-KHEL-" + gisaid_str + "/" + year

    def get_location(self, row):
        return "North America / USA / " + str(row["state"]) if str(row["state"]) != "unknown" else self.default_state

    def get_comment(self, row):
        if str(row['hsn']) in self.priority_lst:
            return self.priority_comment
        else:
            return self.default_comment


def get_hsn(row):
    hsn = str(row["HSN"])
    if len(hsn) == 7:
        return hsn
    else:
        return hsn[:-2]


def get_collection_date(row):
    return str(row["doc"])



def get_sex(row):
    return str(row["sex"]).lower()


def get_priority(row, lst):
    if str(row['hsn']) in lst:
        return 1
    else:
        return 0


