from ms_sql_handler import ms_sql_handler
import pandas as pd
import cx_Oracle as co
import reader
import datetime
from other import add_cols
import json
import os

class demographics_import():

    def __init__(self,cache_path) : #0
        #here need to import json file
        #and used that to store
        demo_cahce= reader.read_json(cache_path+"/data/demographics.json")
        for item in [*demo_cahce] :
            setattr(self,item, demo_cahce[item])
        #and import metric data needed
        self.df_hsn = pd.read_json(cache_path+"/data/sample_metrics.json")
        #df2['year']=df2['year'].astype(int)
        for df_column in self.df_hsn.columns:
            self.df_hsn[df_column]=self.df_hsn[df_column].astype(int)

        
        
    
    def get_lims_demographics(self,hsn,report_dir): #1
        #HSN_WGSRUNDATE_INDEX_
        date= hsn[0].split("_")[1]
        self.wgs_run_date = date
        #base_dict={"FluB_detected":"","FluA_detected":"","FluAH3":"","FluA2009_H1N1":""}
        base_dict={"Influenza B":"","Influenza A":"","Influenza A/H3":"","Influenza 2009 A/H1N1":""}
        #making sure only HSN and cutting other stuff in the array
        hsn = [i.split("_")[0] for i in hsn]
   
        #self.df_hsn = pd.DataFrame(hsn,columns=["hsn"]) #<-could be defined else where, 
        #this will need be a variable i m thinking a jason file that will also feed into msql class
        conn = co.connect(self.lims_connection)

        query="select * from wgsdemographics where HSN in ("+",".join(hsn)+")"
        query_pcr="select HSN,NAME,AMOUNT from FLUWGSDEMO where HSN in ("+",".join(hsn)+") and AMOUNT <> 1000"
        
        self.lims_df = pd.read_sql(query,conn)
        for flu in base_dict:
            self.lims_df[flu] = "Null"


        if not os.path.exists(report_dir+"/"+date):
            os.mkdir(report_dir+"/"+date)

        f= open(report_dir+"/"+date+"/"+date+"_demo.csv","w+")    
        #f.write(",".join(["Accession_Number","Patient_Name","Patient_DOB","Patient_Gender","Ordering_Facility","Specimen_Collection_Date","Race"])+"\n")
        #create txt file of demographical information, for final report
        for sample_df in self.lims_df.values.astype(str).tolist() :
            #5,8,9,13,3,10,
            f.write(",".join([sample_df[0],sample_df[5],sample_df[8],sample_df[9],sample_df[13],sample_df[3],sample_df[10]])+"\n")
        f.close()
            
        
        #print("current lims output after adding new columns")
        #print(self.lims_df.to_string())

        pcr_data_df =pd.read_sql(query_pcr,conn)

        #self.pcr_data={}
        for sample in hsn:
            #self.pcr_data[sample]= base_dict.copy()
            r= pcr_data_df.query("HSN == "+sample).to_dict('records')
            for record in r :
                row_index = self.lims_df[self.lims_df['HSN'] == int(sample)].index.values.astype(int)[0]
           
                self.lims_df.at[row_index,record['NAME']]= record['AMOUNT']
                #self.pcr_data[sample][record['NAME']]=record['AMOUNT']
        

        conn.close()
    
    def format_lims_df(self): #2
        # manipulate sql database to format accepted by the master EXCEL worksheet
        #self.log.write_log("format_lims_DF","Manipulating demographics to database format")
           
        self.lims_df = self.lims_df.rename(columns = self.demo_names)
        self.lims_df["hsn"] = self.lims_df.apply(lambda row: str(row["hsn"]), axis=1)
        #self.log.write_log("format_lims_DF","Done!")


    def merge_dfs(self): #3
        #self.log.write_log("merge_dfs","Merging dataframes")
        self.lims_df['hsn']=self.lims_df['hsn'].astype(int)
      
        self.df = pd.merge(self.lims_df, self.df_hsn, how="inner", on="hsn")
 
       #              self.log.write_log("merge_dfs","Done")
    
    def format_dfs(self): #3 

        #self.log.write_log("format_dfs","Starting")
        # format columns, insert necessary values
        #self.log.write_log("format_dfs","Adding/Formatting/Sorting columns")

        self.df = add_cols(obj=self, \
            df=self.df, \
            col_lst=self.add_col_lst, \
            col_func_map=self.col_func_map)

        # sort/remove columns to match list
        self.df = self.df[self.sample_data_col_order]

        #print(self.df.to_string())
        #self.log.write_log("format_dfs","Done")
    
    def database_push(self): #4
        #self.log.write_log("database_push","Starting")
        self.setup_db()
        df_demo_lst = self.df.values.astype(str).tolist()
        
        self.write_query_tbl1 = (" ").join(self.write_query_tbl1)
   
        self.db_handler.lst_ptr_push(df_lst=df_demo_lst, query=self.write_query_tbl1)

        #need to write this txt file in the results
        #Patient_Last_Name","Patient_First_Name","Patient_DOB","Patient_Gender","Ordering_Facility","Specimen_Collection_Date","Specimen_Source", "Test_Date"
        #self.log.write_log("database_push","Done!`")

    def setup_db(self):
        self.db_handler = ms_sql_handler(self)
        self.db_handler.establish_db()
    



        