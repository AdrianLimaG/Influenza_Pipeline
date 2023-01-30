import os
import subprocess
from reader import get_pandas, read_json
#from Scripts import reader
from other import add_cols, merge_dataframes
from ms_sql_handler import ms_sql_handler


def run_nextclade(sample_name, path_to_results, resource_path, nexclade_output):

   positive_hits={}

   #nextclade_path = path_to_irma+"/resources/nextclade_data/"+what ever flue
   #flu_h1n1pdm_ha,flu_h3n2_ha
   if not os.path.exists(nexclade_output):
    subprocess.run("mkdir "+nexclade_output,shell=True)
  
   resource_path+="/resources"

   for sample in sample_name:
    #run_nextclade
    #temp_hsn = sample.split("_")[0]
    positive_hits[sample]=[]
       
    path_to_sample=path_to_results+"/"+sample+"/amended_consensus"
    #this is where i put the check
    combined = concat_fasta(sample,path_to_sample)
       
    if combined != "FAILED TO ASSEMBLED" :

        subprocess.run("mkdir "+nexclade_output+"/"+sample,shell=True)

        for virus_dataset in ["flu_h3n2_ha", "flu_h1n1pdm_ha"]: #,"flu_yam_ha",,"flu_vic_ha"] : removing this yam_ha isnt working

            nextclade_cmd= resource_path+"/nextclade/nextclade run --in-order --input-dataset "+resource_path+"/nextclade_data/"+virus_dataset+" --output-tsv "+nexclade_output+"/"+sample+"/"+sample+"_"+virus_dataset+".tsv --output-basename "+sample+" "+combined+" --output-all "+nexclade_output+"/"+sample
            #print(nextclade_cmd)
                
            result=subprocess.run(nextclade_cmd, capture_output=True, text=True, shell=True)

            line_n=subprocess.run("wc -l <"+nexclade_output+"/"+sample+"/"+sample+"_"+virus_dataset+".tsv", capture_output=True, text=True, shell=True)
                
            #print(result.stderr) #will need to work this to see which one produce 
            ##nexclade_output+"/"+sample+"/"+sample+"_"+virus_dataset+".tsv 
                    
            if int(line_n.stdout.strip())> 1:
                positive_hits[sample].append(virus_dataset)
            
    else: #failed to have a consensus sequence, need to ask if they want these in the DB
        positive_hits[sample].append("FAILED_TO_ASSEMBLE")
        

    #need to return a dic of dic
    #{samplename:[virusname,virusname too]}else:
   return positive_hits


def concat_fasta(sample_n,path_to_files):

    fastq = os.listdir(path_to_files)
    fastq_s = (" "+path_to_files+"/").join(['',*fastq])
    if os.path.exists(path_to_files):
        if not os.path.exists(path_to_files+"/"+sample_n+"_combined.fasta") :
            output= subprocess.run("cat "+fastq_s+" > "+path_to_files+"/"+sample_n+"_combined.fasta",shell=True)

        return path_to_files+"/"+sample_n+"_combined.fasta"
    else:
        return "FAILED TO ASSEMBLED"


def combine_nextclade_output(nextclade_output_path,samples): 

    #need to loop through all csv files skip frist line
    combined_nextclade = open(nextclade_output_path+"/combined_nextclade.tsv","w+")
    header="seqName\tclade\tqc.overallScore\tqc.overallStatus\ttotalSubstitutions\ttotalDeletions\ttotalInsertions\ttotalFrameShifts\ttotalAminoacidSubstitutions\ttotalAminoacidDeletions\ttotalAminoacidInsertions\ttotalMissing\ttotalNonACGTNs\ttotalPcrPrimerChanges\tsubstitutions\tdeletions\tinsertions\tprivateNucMutations.reversionSubstitutions\tprivateNucMutations.labeledSubstitutions\tprivateNucMutations.unlabeledSubstitutions\tprivateNucMutations.totalReversionSubstitutions\tprivateNucMutations.totalLabeledSubstitutions\tprivateNucMutations.totalUnlabeledSubstitutions\tprivateNucMutations.totalPrivateSubstitutions\tframeShifts\taaSubstitutions\taaDeletions\taaInsertions\tmissing\tnonACGTNs\tpcrPrimerChanges\talignmentScore\talignmentStart\talignmentEnd\tqc.missingData.missingDataThreshold\tqc.missingData.score\tqc.missingData.status\tqc.missingData.totalMissing\tqc.mixedSites.mixedSitesThreshold\tqc.mixedSites.score\tqc.mixedSites.status\tqc.mixedSites.totalMixedSites\tqc.privateMutations.cutoff\tqc.privateMutations.excess\tqc.privateMutations.score\tqc.privateMutations.status\tqc.privateMutations.total\tqc.snpClusters.clusteredSNPs\tqc.snpClusters.score\tqc.snpClusters.status\tqc.snpClusters.totalSNPs\tqc.frameShifts.frameShifts\tqc.frameShifts.totalFrameShifts\tqc.frameShifts.frameShiftsIgnored\tqc.frameShifts.totalFrameShiftsIgnored\tqc.frameShifts.score\tqc.frameShifts.status\tqc.stopCodons.stopCodons\tqc.stopCodons.totalStopCodons\tqc.stopCodons.score\tqc.stopCodons.status\tisReverseComplement\terrors\n"
    
    nextclade_text= []
    temp_clade=""
    nextclade_temp_text =""

    for sample in [*samples]:
        if samples[sample] == ["FAILED_TO_ASSEMBLE"] :
            #means no nextcalde results
            nextclade_text.append(sample+"\tFAILED\tFAILED\tFAILED\n")

        else:    
            i=0
            while i < len(samples[sample]):

                sample_nextcalde = open(nextclade_output_path+"/"+sample+"/"+sample+"_"+samples[sample][i]+".tsv","r")
                lines = sample_nextcalde.readlines()[1]
                temp_l = lines.split("\t")
                #HSN_WGSRUNDATE_INDEX(POS)_PROTIENSEQNUM
                #adding MNUM_RUNNUMBER ---- adding this to thingy
                
                temp_l[0]=temp_l[0]+"_0_0"
                lines = "\t".join(temp_l)

                if i> 0:
                    temp_clade = lines.split("\t")[1] #grabing only clade data
                    
                    last_nextclade = nextclade_text[-1].split("\t")
                    last_nextclade[1]+=";"+temp_clade

                    nextclade_text[-1] = "\t".join(last_nextclade)

                else:
                    nextclade_text.append(lines)

                i+=1

    sample_nextcalde.close()
    combined_nextclade.write(header)
    combined_nextclade.writelines("".join(nextclade_text)) 
    combined_nextclade.close()

    return nextclade_output_path+"/combined_nextclade.tsv"




class nextclade_data_obj():

    def __init__(self,cache_path):

        demo_cahce= read_json(cache_path+"/data/nextclade.json")
        for item in [*demo_cahce] :
            setattr(self,item, demo_cahce[item])
        self.m_num = '1'
        

    def get_nextclade_dfs(self,nc_path):
            # open nextclade path --> pandas dataframe
            #self.log.write_log("get_nextcalde_dfs","Starting function")

            splt = nc_path.split("/")
            parent_folder = splt[-2]
            data = parent_folder.split(".")

            df = get_pandas(nc_path, "WF_4", "nextclade", '\t')
            df = df.rename(columns=self.rename_nc_cols_lst)
            #machine num no longer needed
            #position wil be captured by hsn 
            # add columns
            df = add_cols(obj=self,
                df=df,
                col_lst=self.add_col_lst,
                col_func_map=self.col_func_map)
            
            df.rename(columns={"day_run_num_var":"day_run_num",
                "wgs_run_date_var":"wgs_run_date",
                "machine_num_var":"machine_num"}, inplace=True)
            
            df.fillna("", inplace=True)
            self.df_qc = df[self.nc_qc_cols_lst]
            self.df_results = df[self.nc_results_cols_lst]
            #print(self.df_results)
            #print("------------------------------------")
            #print(self.df_qc)
            #self.log.write_log("get_nextcalde_dfs","complete")
    
    def setup_db(self):
        self.db_handler = ms_sql_handler(self)
        self.db_handler.establish_db()

    def database_push(self):
        #self.log.write_log("database_push","Attemping to connect to db")
        # attempt to connect to database
        self.setup_db()
        #self.log.write_log("database_push","connect to db successful")
        df_qc_update_lst = self.df_qc.values.astype(str).tolist()
        #self.log.write_log("database_push","Pushing information to Run Stats table")
        self.write_query_tbl2 = " ".join(self.write_query_tbl2)
        self.db_handler.lst_ptr_push(df_lst=df_qc_update_lst, query=self.write_query_tbl2)
        self.read_query_tbl2 = " ".join(self.read_query_tbl2)
        all_time_df_qc = self.db_handler.sub_read(query=self.read_query_tbl2)
        
        df_results_final = merge_dataframes(\
            df1=all_time_df_qc, \
            df2=self.df_results, \
            df1_drop=['ID_Table_2', 'percent_cvg', 'avg_depth', 'total_ns'], \
            df_final_drop=['wgs_run_date_y', 'machine_num_y', 'position_y', 'day_run_num_y','wgs_run_date_x', 'machine_num_x', 'position_x', 'day_run_num_x'], \
            join_lst=["hsn", "wgs_run_date", "machine_num", "position", "day_run_num"], \
            join_type='inner')
        

        df_results_final_lst = df_results_final.values.astype(str).tolist()

        if len(df_results_final_lst) == 0:
         #   self.log.write_warning("database_push","Nextclade data from this run has likely already been pushed to the database!")
            raise ValueError("\n-------------------------------------------------------------------------------------------------------------------\
                \nNextclade data from this run has likely already been pushed to the database!\
                \n-------------------------------------------------------------------------------------------------------------------")
       # self.log.write_log("database_push","Updating rows in the results table")
        self.write_query_tbl1 = " ".join(self.write_query_tbl1)
        self.db_handler.lst_ptr_push(df_lst=df_results_final_lst, query=self.write_query_tbl1)
        #self.log.write_log("database_push","Completed")