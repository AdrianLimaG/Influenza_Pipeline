

from merge_seq_data.merge_fastq import merge_seq_fastq
from irma.run_irma import irma_runner
from nextclade.nextcalde import nextclade_runner
from import_demo.import_demo import run_import_demo
import os
import sys

def pipeline(minion_path,sample_sheet_p): #variables, analysis_working_dir, final_out_dir, nextclade_output
    
    res_dir= "/home/ssh_user/FLU_WGS_Sequencing/IRMA/" #this will need a permant address
    res_dir= "/home/ks_khel/Desktop/RES/"
    run_date = sample_sheet_p.split("/")[-1][:-4]
    #adding date to result path
    res_dir=res_dir+run_date+"/"
    dir_path = "/".join(os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]) #path minus scripts 

    nextclade_output= "/home/ssh_user/FLU_WGS_Sequencing/Nextclade"
    nextclade_output=nextclade_output+"/"+run_date
    #Step 1 merge  fasta files
    fastq_paths_dic = merge_seq_fastq(minion_path,sample_sheet_p)

    print("Merging Completing")

    fastq_paths_dic={"2225102_060722_01" : "","2225196_060722_02" :""}

    #Step 2 run irma for allignment
    irma_runner(fastq_paths_dic,dir_path,res_dir) 
    #in the future function after IRMA will return String with where the files have been moved to
    print("IRMA Completing")

    #fastq_paths_dic={'2225196': '/home/ks_khel/Desktop/FLU_DATA/barcode02/2225196_combined.fastq.gz', '2231833': '/home/ks_khel/Desktop/FLU_DATA/barcode04/2231833_combined.fastq.gz', '2225102': '/home/ks_khel/Desktop/FLU_DATA/barcode01/2225102_combined.fastq.gz', '2229929': '/home/ks_khel/Desktop/FLU_DATA/barcode05/2229929_combined.fastq.gz', '2225184': '/home/ks_khel/Desktop/FLU_DATA/barcode03/2225184_combined.fastq.gz'}
    fastq_paths_dic={"2225102_060722_01" : "","2225196_060722_02" :"","2231833_060722_04" :" ","2225184_060722_03" :"",  "2229929_060722_05" :""}
    #Now will need to import demographics from horizon to our local db
    #Step 3 Import Demographics
    #has to be done at this step because this after HSN has been mapped
    run_import_demo(dir_path,[*fastq_paths_dic])


    #Step 4 run nextclade and return hits
    results=nextclade_runner([*fastq_paths_dic],res_dir,dir_path,nextclade_output)
    
    print("Nextclade Completed")


    #next thing to do is to create GISAID Reports 
    #and variant reports
    





if __name__ == "__main__":
    pass
    sample_sheet_p ="/home/ssh_user/FLU_WGS_Sequencing/sample_sheet/060722.csv" #will need a way to pass this differently
    #/home/ks_khel/Desktop/FLU_DATA

    print(sys.argv[0])
    if sys.argv[0] == "" or sys.argv[0] == "/home/ks_khel/Documents/GitHub/Infulenza_Pipeline/Scripts/flu_pipeline.py":
        input_path= input("Please enter the path to MinION data:     ")

        pipeline(input_path,sample_sheet_p)
    else:
        #pipeline(sys.argv[0])
        pipeline("/home/ssh_user/FLU_WGS_Sequencing/run_data/060722",sample_sheet_p)


#HSN****,PLATEPOS,INDEX**

