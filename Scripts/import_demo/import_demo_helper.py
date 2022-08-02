from ms_sql_handler import ms_sql_handler
import pandas as pd
import cx_Oracle as co
import reader
import datetime
from other import add_cols

class demographics_import():

    def __init__(self,cache_path) : #0
        #here need to import json file
        #and used that to store
        demo_cahce= reader.read_json(cache_path+"/data/demographics.json")
        for item in [*demo_cahce] :
            setattr(self,item, demo_cahce[item])
        
    
    def get_lims_demographics(self,hsn): #1
        #HSN_WGSRUNDATE_INDEX_
        #making sure only HSN and cutting other stuff in the array
        hsn = [i.split("_")[0] for i in hsn]    
        self.df_hsn = pd.DataFrame(hsn,columns=["hsn"])
        #this will need be a variable i m thinking a jason file that will also feed into msql class
        conn = co.connect(self.lims_connection)

        query="select * from wgsdemographics where HSN in ("+",".join(hsn)+")"

        self.lims_df = pd.read_sql(query,conn)

        conn.close()
    
    def format_lims_df(self): #2
        # manipulate sql database to format accepted by the master EXCEL worksheet
        #self.log.write_log("format_lims_DF","Manipulating demographics to database format")
        self.lims_df = self.lims_df.rename(columns = self.demo_names)
        self.lims_df["hsn"] = self.lims_df.apply(lambda row: str(row["hsn"]), axis=1)
        #self.log.write_log("format_lims_DF","Done!")


    def merge_dfs(self): #3
        #self.log.write_log("merge_dfs","Merging dataframes")
        self.df = pd.merge(self.lims_df, self.df_hsn, how="right", on="hsn")
        #self.log.write_log("merge_dfs","Done")
    
    def format_dfs(self, runId): #3 

        #self.log.write_log("format_dfs","Starting")
        # get the date for wgs_run_date column
        #not run ID but possibly file name
        self.wgs_run_date = datetime.datetime.strptime(runId, '%Y-%m-%d').strftime("%m/%d/%Y")

        # format columns, insert necessary values
        #self.log.write_log("format_dfs","Adding/Formatting/Sorting columns")

        self.df = add_cols(obj=self, \
            df=self.df, \
            col_lst=self.add_col_lst, \
            col_func_map=self.col_func_map)

        # sort/remove columns to match list
        self.df = self.df[self.sample_data_col_order]
        #self.log.write_log("format_dfs","Done")
    
    def database_push(self): #4
        #self.log.write_log("database_push","Starting")
        self.setup_db()
        df_demo_lst = self.df.values.astype(str).tolist()
        df_table_col_query = "(" + ", ".join(self.df.columns.astype(str).tolist()) + ")"
        self.write_query_tbl1 = self.write_query_tbl1.replace("{df_table_col_query}", df_table_col_query)
        self.db_handler.lst_ptr_push(df_lst=df_demo_lst, query=self.write_query_tbl1)
        #self.log.write_log("database_push","Done!`")

    def setup_db(self):
        self.db_handler = ms_sql_handler(self)
        self.db_handler.establish_db()



        