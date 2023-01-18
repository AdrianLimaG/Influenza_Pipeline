

from WF_0_merge_seq_data.merge_fastq import merge_seq_fastq
from WF_1_irma.run_irma import irma_runner
from WF_3_nextclade.nextcalde import nextclade_runner
from WF_2_import_demo.import_demo import run_import_demo
from WF_5_final_report.WF_5_final_report import create_final_report, move_fasta_files, create_alignment_file, create_phylogentic_tree, clean_run_files
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
        print("Merging Completing")
        
        # #Step 1 run irma for allignment
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
    
        #create Final Report for Epis
        create_final_report(run_date,self.nextclade_output,results,self.final_results_dir)

        #Move fasta files into result dir
        move_fasta_files([*fastq_paths_dic],self.res_dir,self.final_results_dir, run_date )

        #create alignment files per flu type tested
        create_alignment_file(self.dir_path,self.res_dir,results,self.final_results_dir,run_date)
        print("Alignment Files Created")
        #create tree based on alignment files
        create_phylogentic_tree(self.final_results_dir,run_date)
        print("Phlyogenetic Tree Built")

        #delete barcodes + IRMA data + Nextclade
        clean_run_files(minion_path,self.res_dir,run_date,self.nextclade_output)


        print("PipeLine Finished!")



if __name__ == "__main__":
    
    dir_path = "/".join(os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]) #path minus scripts 
    
    #print(sys.argv)
    input_path = sys.argv[1]
    sample_sheet_p = sys.argv[2]
    print(input_path)
    print("-----------")
    print(sample_sheet_p)
    print("-----------")
 
    influenza_pipeline = flu_pipeline(dir_path)
       
    influenza_pipeline.run_flu_pipeline(input_path,sample_sheet_p)







