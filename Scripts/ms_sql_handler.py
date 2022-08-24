from sqlalchemy import create_engine
import re
import time
import sys
import pandas as pd
import re



#Table = "dbo.Results"
#Table = "dbo.Run_Stats"

class ms_sql_handler():
    # constructor:
    def __init__(self, obj):
        self.sql_user = obj.sql_user
        self.sql_pass = obj.sql_pass
        self.sql_server = obj.sql_server
        self.sql_db = obj.sql_db
        #self.logger = obj.logger
       
        self.avg_depth_cutoff = obj.avg_depth_cutoff
        self.percent_cvg_cutoff = obj.percent_cvg_cutoff
        self.engine = None
#        self.log = Script_Logger("SQL_Handler")
#        self.log.start_log("Starting")

    # methods
    # Create
    def establish_db(self):
        try:
            self.engine = create_engine('mssql+pyodbc://' + self.sql_user + ':' + self.sql_pass + '@' + self.sql_server + '/' + self.sql_db + '?driver=ODBC+Driver+17+for+SQL+Server')
        except Exception as e:
            print(e)
            #self.logger.critical(self.log_name + ": Issue in connection to mssql database")
            time.sleep(20)
            sys.exit

    # Clear (DELETE)
    def clear_db(self):
        with self.engine.connect() as conn:
            #self.logger.info("refresh: deleting all information from database")
            query = "DELETE FROM dbo.Results"
            res = conn.execute(query)
            query = "DELETE FROM dbo.Run_Stats"
            res = conn.execute(query)
        #self.logger.info("refresh: database_clear finished!")
    

    # Read
    def ss_read(self, query=None):
        with self.engine.connect() as conn:
            df = pd.read_sql(query, con=self.engine)
        return df

    def sub_read(self, query=None):
        local_query = query.replace("{avg_depth_cutoff}", str(self.avg_depth_cutoff))
        local_query = local_query.replace("{percent_cvg_cutoff}", str(self.percent_cvg_cutoff))
        df = pd.read_sql(local_query, con=self.engine)
        return df

    def sub_lst_read(self, query=None, lst=None):
        hsn_query = "(" + ", ".join(lst) + ")"
        local_query = query.replace("{hsn_query}", hsn_query)
        # create dataframe that has structure of table_1 but only includes HSNs from GISAID.xlsx
        df = pd.read_sql(local_query, con=self.engine)
        return df


    # Write
    def to_sql_push(self, df=None, tbl_name=None, u_if_exists='append', u_index=False):
        df.to_sql(tbl_name, self.engine, if_exists=u_if_exists, index=u_index)

    def lst_push(self, df_lst=None, df_cols=None):
        """
        push results to database with the replace strategy
        @params:
            df_lst      - Required  : nested list that holds values to fill (List[List[]])
            df_cols     - Required  : list of column names for supplied list (List[])
        """
        with self.engine.connect() as conn:
            for i in range(len(df_lst)):
                df_lst[i] = format_lst(df_lst[i])
                df_lst_query = "(" + ", ".join(df_lst[i]) + ")"
                query = f"""INSERT INTO dbo.Run_Stats {df_cols} VALUES {df_lst_query}"""
                res = conn.execute(query)

    def lst_ptr_push(self, df_lst=None, query=None, full=False, df=None):
        """
        push results to database with the replace (list/pointer) strategy
        @params:
            df_lst      - Required  : nested list that holds values to fill ([[x],[y]])
            query       - Required  : stored in cache, generic string with placeholders (Str)
            full        - Optional  : True if full dataframe being pushed. For use with refresh or ouside_lab scripts (Bool)
            df          - Optional  : dataframe of values to clean. For use with refresh or outside_lab (Pandas DataFrame)
        """
        local_query = query
        # connect to db
        print("trying to connect to sql")
        with self.engine.connect() as conn:
            print("connection passed!")
            # generate a distinct query for every row, where query stores the
            # generic value
            for i in range(len(df_lst)):
                # if outside_lab or refresh, we are using full excel file, replace
                # as needed
                if full:
                    df_table_col_lst = list(df.columns)
                    # remove any columns/entries that are 'nan' or 'None'
                    element = 0
                    df_ctr = 0
                    while element < len(df_table_col_lst):
                        x = str(df.iloc[i, df_ctr])
                        if x == "nan" or x == "None" or x == "extraction only, WGS":
                            self.log.write_warning("removing None/nan","This is being removed   "+x)
                            del df_table_col_lst[element]
                            df_ctr += 1
                        else:
                            self.log.write_log("This is being kept","Kepping this  "+x)
                            element += 1
                            df_ctr += 1
                    # now, create the col list
                    df_table_col_query = "(" + ", ".join(df_table_col_lst) + ") "
                
                # new_query stores distinct query for the corresponding row
                new_query = local_query

                if full:
                    new_query = new_query.replace("{df_table_col_query}", df_table_col_query)
                # find all unique occurrances of '{/d}' and add them to a list
                query_track = list(set(re.findall("({.*?})", new_query)))
                
                # replace all occurrances of '{\d}' with the corresponding
                # values in the df_lst
                try:
                    for item in query_track:

                        new_query = new_query.replace(item, df_lst[i][int(item[1:-1])])

                except IndexError:
                    pass
                # if outside_lab or refresh, we are using full excel file, replace
                # missing data as needed
                if full:
                    new_query = new_query.replace(", CAST(\'nan\' AS DATE)", "")
                    new_query = new_query.replace(", CAST(\'None\' AS DATE)", "")
                    new_query = new_query.replace(", \'nan\'", "")
                    new_query = new_query.replace(", \'None\'", "")
                    new_query = new_query.replace(", nan", "")
                    new_query = new_query.replace(", None", "")
                    new_query = new_query.replace(", \'extraction only, WGS\'", "")
                    new_query = new_query.replace("other", "OT")
                    new_query = new_query.replace("(, ", "(")
                    new_query = new_query.replace(" KS", " Kansas")
                new_query = new_query.replace("= ,", "= NULL,")
                new_query = new_query.replace("= '',", "= NULL,")
                new_query = new_query.replace('= "",', '= NULL,')
                new_query = new_query.replace("= \'None\',", "= NULL,")
                new_query = new_query.replace("CAST('nan' AS DATE)", "NULL")
                new_query = new_query.replace("luke's", 'lukes')
                new_query = new_query.replace("'None'", "NULL")
                res = conn.execute(new_query)


def format_lst(lst):
    #iterate through the list
    for i in range(len(lst)):
        try:

            # we only want to add quotes to entries that are not null
            if not re.search(lst[i], "null"):
                # attempt to convert them to float
                try_float = float(lst[i])
        except ValueError:
            # valueError raised if float convert is attempted with string
            lst[i] = f"'{lst[i]}'"
    return lst