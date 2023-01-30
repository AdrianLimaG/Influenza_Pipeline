

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
        #fastq_paths_dic = merge_seq_fastq(minion_path,sample_sheet_p)
        print("Merging Completing")
        
        # #Step 1 run irma for allignment
        #irma_runner(fastq_paths_dic,self.dir_path,self.res_dir) 
        print("IRMA Completing")
        fastq_paths_dic ={'2405320_012723_81':'',
                            '2405321_012723_82':'',
                            '2403806_012723_65':'',
                            '2403813_012723_68':'',
                            '2404051_012723_66':'',
                            '2404061_012723_67':'',
                            '2404987_012723_69':'',
                            '2405011_012723_78':'',
                            '2405013_012723_79':'',
                            '2405020_012723_72':'',
                            '2405312_012723_73':'',
                            '2405313_012723_74':'',
                            '2405314_012723_75':'',
                            '2405315_012723_76':'',
                            '2405316_012723_77':'',
                            '2405317_012723_80':'',
                            '2405314_012723_75':'',
                            '2405315_012723_76':'',
                            '2405316_012723_77':'',
                            '2405317_012723_80':'',
                            '2405322_012723_83':'',
                            '2405324_012723_84':'',
                            '2411612_012723_89':'',
                            '2411613_012723_90':'',
                            '2411615_012723_91':'',
                            '2411622_012723_92':'',
                            '2411623_012723_93':'',
                            '2411901_012723_85':'',
                            '2411903_012723_86':'',
                            '2411904_012723_87':'',
                            '2411906_012723_88':'',
                            '2411907_012723_94':'',
                            '2417077_012723_95':'',
                            '2419133_012723_96':'',
                            '2411906_012723_88':'',
                            '2411907_012723_94':'',
                            '2417077_012723_95':'',
                            '2419133_012723_96':'',
                            '2419136_012723_70':'',
                            '2419137_012723_71':''}
        #Step 2 Import Demographics
        #run_import_demo(self.dir_path,[*fastq_paths_dic],self.final_results_dir)

        #Step 3 run nextclade and return hits
        #results=nextclade_runner([*fastq_paths_dic],self.res_dir,self.dir_path,self.nextclade_output)
        print("Nextclade Completed")
        results ={'2405320_012723_81': ['flu_h3n2_ha'], '2405321_012723_82': ['flu_h3n2_ha'], '2403806_012723_65': ['flu_h3n2_ha'], '2403813_012723_68': ['flu_h3n2_ha'], '2404051_012723_66': ['flu_h3n2_ha'], '2404061_012723_67': ['flu_h3n2_ha'], '2404987_012723_69': ['flu_h3n2_ha'], '2405011_012723_78': ['flu_h3n2_ha'], '2405013_012723_79': ['flu_h3n2_ha'], '2405020_012723_72': ['flu_h3n2_ha'], '2405312_012723_73': ['flu_h3n2_ha'], '2405313_012723_74': ['flu_h3n2_ha'], '2405314_012723_75': ['flu_h3n2_ha'], '2405315_012723_76': ['flu_h3n2_ha'], '2405316_012723_77': ['flu_h3n2_ha'], '2405317_012723_80': ['flu_h3n2_ha'], '2405322_012723_83': [], '2405324_012723_84': ['flu_h1n1pdm_ha'], '2411612_012723_89': ['flu_h3n2_ha'], '2411613_012723_90': ['flu_h3n2_ha'], '2411615_012723_91': ['flu_h3n2_ha'], '2411622_012723_92': ['flu_h3n2_ha'], '2411623_012723_93': ['flu_h1n1pdm_ha'], '2411901_012723_85': ['flu_h3n2_ha'], '2411903_012723_86': ['flu_h3n2_ha'], '2411904_012723_87': ['flu_h3n2_ha'], '2411906_012723_88': ['flu_h3n2_ha'], '2411907_012723_94': ['flu_h3n2_ha'], '2417077_012723_95': ['flu_h1n1pdm_ha'], '2419133_012723_96': ['flu_h1n1pdm_ha'], '2419136_012723_70': ['flu_h1n1pdm_ha'], '2419137_012723_71': ['flu_h3n2_ha']}
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
    sample_sheet_p ="/FLU_WGS_Sequencing/sample_sheet/060722.csv" #will need a way to pass this differently
 
    print(sys.argv)
    input_path = sys.argv[1]
    sample_sheet_p = sys.argv[2]
    print(input_path)
    print("-----------")
    print(sample_sheet_p)
    print("-----------")
 
    influenza_pipeline = flu_pipeline(dir_path)
       
    influenza_pipeline.run_flu_pipeline(input_path,sample_sheet_p)







