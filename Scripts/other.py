import datetime
import pandas as pd
import numpy as np
import os
import re

def merge_dataframes(df1=None, df2=None, df1_drop=None, df_final_drop=None, join_lst=None, join_type=None):
    # df1 from qc table, may have duplicate hsns.  Remove common columns between
    # the two dataframes. df2 from results table
    join_dict = {}
    for col in join_lst:
        join_dict[col] = 'str'
    df1 = df1.astype(join_dict)
    df2 = df2.astype(join_dict)

    df1.drop(labels=df1_drop, axis=1, inplace=True)

    #df1.drop_duplicates(subset=['hsn'], inplace=True)
    df_final = df2.merge(df1, how=join_type, on='hsn')
    df_final.drop(labels=df_final_drop, axis=1, inplace=True)

    return df_final


def add_cols(obj=None, df=None, col_lst=None, col_func_map=None):
    # iterate through all columns that should be in final df
    for k in col_lst:
        # if the column appears in the function mapping,
        if k in col_func_map.keys():
            # get the pointer to the func/value associated with the column
            v = col_func_map[k]
            try:
                # try to get additional value to run apply function with
                val = v[1]
                try:
                    val = getattr(obj, v[1])

                    # try to catch v[1] as an object variable
                    df[k] = df.apply(lambda row: globals()[v[0]](row, val), axis=1)
                except Exception:
                    # use v[1] as a constant argument to the function
                   
                    df[k] = df.apply(lambda row: globals()[v[0]](row, v[1]), axis=1)
            # no additional variables to supply to apply function
            except IndexError:
                # try using the value as a function
                try:
                    df[k] = df.apply(lambda row: globals()[v[0]](row), axis=1)
                # try using the value as a variable
                except Exception:
                    val = getattr(obj, v[0])
                    df[k] = val
        # if column not in mapping, insert empty column with appropriate
        # name into the dataframe
        else:
            # try to insert the column
            try:
                df.insert(0, k, None)
            # ValueError raised if column already exists in dataframe
            except ValueError:
                pass

    return df





def format_date(row, colName):
    if (isinstance(row[colName], pd.Timestamp)) or\
        (isinstance(row[colName], datetime.datetime)):
        if (not pd.isna(row[colName])):
            return row[colName].strftime("%m/%d/%Y")
        else:
            return np.nan
    else:
        return np.nan
        

def get_today(row):
    return datetime.datetime.today().strftime("%Y-%m-%d")

def get_pos(id):
    pos_dict = {"A":1, "B":2, "C":3, "D":4, "E":5, "F":6, "G":7, "H":8}
    pos = (int(id[-1])*8 - 8) + pos_dict[id[0]]
    return pos


def parse_seq_id(row, arg):
    #print(row)
    #print(arg)
    try:
        seq_id = str(row["Sequence name"]).split("_")
    except:
        seq_id = str(row['seqName']).split("_")
    # if split didn't find matches, it is dealing with folder, should
    # be split by ".", also has different indexes for values
    if len(seq_id) == 1: #this won't happend anymore
        # WF 3
        seq_id = str(row["seqName"]).split(".")
        if arg == "hsn":
            return int(seq_id[0][0:7])
        elif arg == "m_num":
            return int(seq_id[1][-2:])
        elif arg == "pos":
            return int(seq_id[4][-2:])
        elif arg == "run_num":
            return int(seq_id[3])
        elif arg == "date":
            return seq_id[2]
        else:
            raise ValueError("Bad argument to parse_seq_id --> folder")
    else:
        # WF_4, WF_5
        #HSN_WGSRUNDATE_INDEX(POS)_PROTIENSEQNUM_MNUM_RUNNUM
        if arg == "hsn":
            return int(seq_id[0])
        elif arg == "m_num":
            return int(seq_id[4])
        elif arg == "pos":
            return int(seq_id[2])
        elif arg == "run_num":
            return int(seq_id[5])
        elif arg == "date":
            return seq_id[1]
        else:
            raise ValueError("Bad argument to parse_seq_id --> file")


def extract_hsn(row):
    hsn = str(row["Sample ID"])
    if len(hsn) == 7:
        return hsn
    return hsn[:-2]


def format_str_cols(df):
    df.columns = [str(col) for col in list(df.columns)]
    return df


def format_sex(row, ber=False):
    col = "sex"
    if ber:
        col = 'Patient_Gender'
    
    if pd.isna(row[col]) or str(row[col]).upper() == "" or str(row[col]).upper() == "UNKNOWN" or str(row[col]).upper() == "U":
        return "Unknown"
    elif str(row[col]).upper() == "M":
        return "Male"
    elif str(row[col]).upper() == "F":
        return "Female"
    else:
        return str(row[col]).capitalize()


def format_facility(row, facility_replace_dict):
    if pd.isna(row['facility']):
        return None
    elif row['facility'] == "":
        return None
    else:
        facility = str(row['facility']).lower()
        for key in facility_replace_dict.keys():
            facility = facility.replace(key, facility_replace_dict[key])
        return facility.lower()


def parse_category(row, parse_category_dict):
    facility = str(row['facility']).lower()
    for key in parse_category_dict.keys():
        if re.search(key, facility):
            return parse_category_dict[key]
    return None


def format_race(row):
    if pd.isna(row['race']) or row['race'] == "U":
        return "Unknown"
    elif row['race'] == "":
        return "Unknown"
    elif str(row['race']).upper() == "W":
        return "White"
    else:
        return str(row['race'])


def format_source(row):
    source = str(row['source']).lower()
    if len(source) > 2:
        if source == "nasopharyngeal":
            return "NP"
        elif source == "sputum/saliva":
            return "SV"
        else:
            return "OT"
    else:
        return row['source']


def add_cols_by_name(df, lst):
    curr_cols = list(df.columns)
    for col in lst:
        if not (col in curr_cols):
            df.insert(0, col, np.nan)
    return df


def format_f_name(row):
    if pd.isna(row['name']):
        return None
    elif row['name'].strip() == "":
        return None
    else:
        full_name = str(row["name"])
        names = full_name.split()
        f_name = names[0].capitalize()
        f_name = f_name.replace("'","''")
        return f_name


def format_l_name(row, lst):
    if pd.isna(row['name']):
        return None
    elif row['name'].strip() == "":
        return None
    else:
        full_name = str(row["name"])
        names = full_name.split()
        for item in lst:
            if item == names[-1].lower():
                return names[-2].capitalize() + ", " + names[-1].upper()
        l_name = names[-1].capitalize()
        l_name = l_name.replace("'", "''")
        return l_name


def drop_cols(df, lst):
    
    curr_cols = list(df.columns)
    for col in curr_cols:
        if not (col in lst):
            df = df.drop(columns = col)
    return df


def get_age(row):
    try:
        born = datetime.datetime.strptime(row["dob"], "%m/%d/%Y").date()
    except Exception:
        born = row["dob"].to_pydatetime().date()

    try:
        tested = row['doc'].to_pydatetime().date()
    except Exception:
        tested = datetime.datetime.strptime(row['doc'], "%m/%d/%Y").date()
    if pd.isnull(born) or pd.isnull(tested):
        return -1
    days_in_year = 365.2425   
    age = int((tested - born).days / days_in_year)
    return age


def format_state(row, state_abbrev):
    if (not pd.isna(row['state'])):
        return state_abbrev[str(row["state"])]
    else:
        return "unknown"


def parse_path(row, path):
    new_path = path + "/" + row["seqName"]
    if not os.path.exists(new_path):
        raise ValueError("The parser generated a path to a fasta file that is not valid!!")
    return new_path

def get_gisaid(row):
    if np.isnan(row["gisaid_num"]) or pd.isna(row["gisaid_num"]):
        return ""
    else:
        return "KS-KHEL-" + str(int(row["gisaid_num"]))


def get_name(row):
    return str(row['f_name']) + " " + str(row['l_name'])


def unkwn(row, col):
    if (str(row[col]) == "nan") or (str(row[col]) == "Nan"):
        return ""
    else:
        return row[col]


def cap_all(row, col):
    splt = str(row[col]).split()
    newlst = []
    for word in splt:
        newlst.append(word.capitalize())
    return " ".join(newlst)


def get_priority(row, lst):
    if row['hsn'] in lst:
        return 1
    else:
        return 0


def check_reportable(row, cutoff):
    if row['neg_pass'] and row['pos_pass'] and row['percent_cvg'] >= cutoff:
        return 1
    else:
        return 0


def replace_shortcut(path):
    base_path = "//kdhe/dfs/LabShared"
    path = path.replace("\\", "/")
    if base_path == path[0:20]:
        return path
    if path[3:26] != "Molecular Genomics Unit":
        # return None
        raise ValueError("The supplied path shouldn't be passed to 'replace_shortcut()'!")
    extension = path[3:]
    path = base_path + "/" + extension
    return path