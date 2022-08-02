import time
import json
import pandas as pd
import sys


def get_pandas(path, log_name, workbook_name, separator):
    try:
        ext = path.split(".")[-1]
        if ext == "csv" or ext == "tsv":
            df = pd.read_csv(path, sep=separator)
        else:
            df = pd.read_excel(path)
    except Exception as o:
        print("\nThere is an issue opening the " + workbook_name + " workbook!")
        print(o)
        time.sleep(10)
        sys.exit()
    print(" Done!\n")
    return df


def read_json(path):
    with open(path, 'r') as json_cache:
        full_cache = json.load(json_cache)
    return full_cache


def read_txt(path):
    with open(path, 'r') as f:
        lst = f.readlines()
    return lst