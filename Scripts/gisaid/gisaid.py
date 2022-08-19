
from gisaid.gisaid_helper import gisaid_obj

def run_gisaid(cache_path):
    print("\n================================\nGISAID Report Script\n================================\n\n")

    # import relevant data from json file
    data_obj = gisaid_obj(cache_path)

    data_obj.scan_db()
    data_obj.get_gisaid_df()
    data_obj.compile_fasta(run_id)
    data_obj.compile_gisaid()
    data_obj.make_fasta_file()
    data_obj.make_gisaid_file(run_id)
    data_obj.database_push()


    print("\n================================\nSUCCESS - END OF SCRIPT\n================================\n\n")
