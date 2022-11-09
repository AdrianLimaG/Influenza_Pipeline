

from WF_0_merge_seq_data.merge_fastq import merge_seq_fastq
from WF_1_irma.run_irma import irma_runner
from WF_3_nextclade.nextcalde import nextclade_runner
from WF_2_import_demo.import_demo import run_import_demo
from WF_5_final_report.WF_5_final_report import create_final_report, move_fasta_files
import os
import sys
import pandas as pd
import reader



class flu_pipeline() :

    def __init__(self,cache_path) :
        self.dir_path= cache_path
        demo_cahce= reader.read_json(cache_path+"/data/pipeline_resources.json")
        for item in [*demo_cahce] :
            setattr(self,item, demo_cahce[item])


    def run_flu_pipeline(self,minion_path,sample_sheet_p): #variables, analysis_working_dir, final_out_dir, nextclade_output
        run_date = sample_sheet_p.split("/")[-1][:-4]
    
        self.res_dir= self.res_dir+run_date
        self.nextclade_output= self.nextclade_output+"/"+run_date

        #Step 0 merge  fasta files
        fastq_paths_dic = merge_seq_fastq(minion_path,sample_sheet_p)
        print("Merging Completing")#,'2208336_060722_31':''
        #fastq_paths_dic={'2200158_060722_18':'','2200185_060722_20':'','2200209_060722_19':'','2200229_060722_17':'','2201125_060722_22':'','2201150_060722_21':'','2202121_060722_23':'','2204939_060722_24':'','2204972_060722_25':'','2204990_060722_26':'','2205005_060722_29':'','2205018_060722_28':'','2205025_060722_27':'','2208361_060722_32':'','2208364_060722_30':'','2225102_060722_01':'','2225184_060722_03':'','2225196_060722_02':'','2229916_060722_06':'','2229929_060722_05':'','2231833_060722_04':'','2231954_060722_07':'','2233275_060722_08':'','2233282_060722_10':'','2234116_060722_09':'','2234922_060722_11':'','2234924_060722_12':'','2239358_060722_13':'','2240979_060722_14':'','2244254_060722_15':'','2245770_060722_16':''}
        #Step 1 run irma for allignment
        irma_runner(fastq_paths_dic,self.dir_path,self.res_dir) 
        print("IRMA Completing")
        
        #Step 2 Import Demographics
        run_import_demo(self.dir_path,[*fastq_paths_dic],self.final_results_dir)

        #Step 3 run nextclade and return hits
        results=nextclade_runner([*fastq_paths_dic],self.res_dir,self.dir_path,self.nextclade_output)
        print("Nextclade Completed")

        #results were supose to be passed to Step 4 to create GISAID upload report

        #Step 5 create final Report and move fast files all together
        #I will create a script to output a tsv of hsn, and clade hit
        #That needs to be moved back to the main network from Analysis PC
        #self.final_results_dir This is were i will dump the value of results into a file and the aligned sequences
    
        create_final_report(run_date,self.nextclade_output,results,self.final_results_dir)

        move_fasta_files([*fastq_paths_dic],self.res_dir,self.final_results_dir, run_date )


        print("PipeLine Finished!")



def pipeline(minion_path,sample_sheet_p): #variables, analysis_working_dir, final_out_dir, nextclade_output
    
    res_dir= "/home/ssh_user/FLU_WGS_Sequencing/IRMA/" #this will need a permant address

    run_date = sample_sheet_p.split("/")[-1][:-4]

    dir_path = "/".join(os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]) #path minus scripts 
    res_dir+=run_date
    nextclade_output= "/home/ssh_user/FLU_WGS_Sequencing/Nextclade"
    nextclade_output=nextclade_output+"/"+run_date

    #Step 1 merge  fasta files
    #fastq_paths_dic = merge_seq_fastq(minion_path,sample_sheet_p)

    print("Merging Completing")

    fastq_paths_dic={"2225102_060722_01" : "","2225196_060722_02" :""}

    #Step 2 run irma for allignment
    #irma_runner(fastq_paths_dic,dir_path,res_dir) 
    #in the future function after IRMA will return String with where the files have been moved to
    print("IRMA Completing")

    #fastq_paths_dic={'2225196': '/home/ks_khel/Desktop/FLU_DATA/barcode02/2225196_combined.fastq.gz', '2231833': '/home/ks_khel/Desktop/FLU_DATA/barcode04/2231833_combined.fastq.gz', '2225102': '/home/ks_khel/Desktop/FLU_DATA/barcode01/2225102_combined.fastq.gz', '2229929': '/home/ks_khel/Desktop/FLU_DATA/barcode05/2229929_combined.fastq.gz', '2225184': '/home/ks_khel/Desktop/FLU_DATA/barcode03/2225184_combined.fastq.gz'}
   # fastq_paths_dic={"2225102_060722_01" : "","2225196_060722_02" :"","2231833_060722_04" :" ","2225184_060722_03" :"",  "2229929_060722_05" :""}
    
    #Now will need to import demographics from horizon to our local db
    #Step 3 Import Demographics
    #has to be done at this step because this after HSN has been mapped
    run_import_demo(dir_path,[*fastq_paths_dic])


    #Step 4 run nextclade and return hits
    #results=nextclade_runner([*fastq_paths_dic],res_dir,dir_path,nextclade_output)
    
    print("Nextclade Completed")

    #next thing to do is to create GISAID Reports #epiFLU seems to be still down
    #and variant reports
    



if __name__ == "__main__":
    
    dir_path = "/".join(os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]) #path minus scripts 
    sample_sheet_p ="/home/ssh_user/FLU_WGS_Sequencing/sample_sheet/060722.csv" #will need a way to pass this differently
    #/home/ks_khel/Desktop/FLU_DATA
    #print(sys.argv)
    input_path = sys.argv[1]
    sample_sheet_p = sys.argv[2]
    print(input_path)
    print("-----------")
    print(sample_sheet_p)
    #print(sys.argv[0])
    influenza_pipeline = flu_pipeline(dir_path)
       
    influenza_pipeline.run_flu_pipeline(input_path,sample_sheet_p)

   # if sys.argv[0] == "" or sys.argv[0] == "/home/ssh_user/Documents/GitHub/Infulenza_Pipeline/Scripts/flu_pipeline.py":
    #    input_path= input("Please enter the path to MinION data:     ")

        #pipeline(input_path,sample_sheet_p)
        
     #   influenza_pipeline = flu_pipeline(dir_path)
        
      #  influenza_pipeline.run_flu_pipeline(input_path,sample_sheet_p)







