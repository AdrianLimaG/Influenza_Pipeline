import os
import subprocess
from datetime import date
import shutil



def create_final_report(run_date,path_to_nextclade,nextclade_hits,result_output):
    
    today = date.today().strftime("%b-%d-%Y")
    hits={}
    #could be wirten better
    for key in [*nextclade_hits] :
        
        hits[key.split("_")[0]] = nextclade_hits[key]


    header=["Accession_Number","Patient_Name","Patient_DOB","Patient_Gender","Ordering_Facility", \
        "Specimen_Collection_Date","Race","Influenza_Type","Lineage_ID","Nextclade_Score","WGS_Run_Date","Report_Generation_Date"]
    #report_dir+"/"+date+"/"+date+"_demo.txt" this is the demo stuff
    #"Accession_Number","Patient_Name","Patient_DOB","Patient_Gender","Ordering_Facility","Specimen_Collection_Date","Race"

    #create report 
    report = open(result_output+"/"+run_date+"/"+run_date+"_report.tsv","w+")
    report.write("\t".join(header)+"\n")

    #read in demo_file
    demo= open(result_output+"/"+run_date+"/"+run_date+"_demo.csv","r")

    #read_nextclade_file
    nextclade_file = open(path_to_nextclade+"/combined_nextclade.tsv")
    nextclade_dict={}
    for l in nextclade_file.readlines():
        l=l.strip()
        l=l.split("\t")
        nextclade_dict[l[0].split("_")[0]]=l[1],l[3]
    print(nextclade_dict)

   
    for line in demo.readlines() :
        
        #"Influenza_Type","Lineage_ID","Nextclade_Score","WGS_Run_Date","Report_Generation_Date"
        line = line.strip().split(",")
        hsn = line[0]
        if len(hits[hsn])>0:
                                    #flu type                  #linageID            #score
            report_line = line + hits[hsn][:] + [nextclade_dict[hsn][0]] + [nextclade_dict[hsn][1]]+ [run_date] + [today]
        else:
            report_line = line + ["NO HITS"] + ["NA"] + ["NA"]+ [run_date] + [today]

        report.write("\t".join(report_line)+"\n")
    demo.close()
    report.close()
    #clean up demo report
    subprocess.run("rm "+result_output+"/"+run_date+"/"+run_date+"_demo.csv",shell=True)
    print("report created")


def move_fasta_files(fastq_samples,path_to_irma,result_output_dir,runD):

    os.mkdir(result_output_dir+"/"+runD+"/fasta_files")

    for sample in fastq_samples:
        temp_hsn = sample.split("_")[0]

        subprocess.run("cp "+path_to_irma+"/"+sample+"/amended_consensus/"+sample+"_combined.fasta "+result_output_dir+"/"+runD+"/fasta_files/"+temp_hsn+".fasta", shell=True)

    
    print("Fasta Files have been copied")


#need two new function to create phylo genetic tree
#one to create alignmnet file
    #copy all XXX_4.fa files into a temp dir then run it
#output dir should include runID
def create_alignment_file (path_to_resources,path_to_irma,nextclade_hsn_file_name,output_dir,runD):
    
    output_dir+="/"+runD
    #path_to_irma+="/"+runD
    for virus_dataset in ["flu_h3n2_ha", "flu_h1n1pdm_ha"]:

        os.mkdir(output_dir+"/temp_alignment_"+virus_dataset)
        os.mkdir(output_dir+"/fasta_alignment_"+virus_dataset)

        #get all files in one place
        for sample in [*nextclade_hsn_file_name]:
            #hsn = sample.split('_')[0]
            if virus_dataset in nextclade_hsn_file_name[sample] :
                subprocess.run("cp "+path_to_irma+"/"+sample+"/amended_consensus/*_4.fa "+output_dir+"/temp_alignment_"+virus_dataset,shell=True)
        
    #nextalign run --input-ref=nextclade-master/data/flu_h3n2_ha/reference.fasta 
    #--genemap=nextclade-master/data/flu_h3n2_ha/genemap.gff --output-all=output_nextali/ temp_align/*.fa

        subprocess.run(". $CONDA_PREFIX/home/ssh_user/miniconda3/etc/profile.d/conda.sh && conda activate nextali && nextalign run --input-ref="+path_to_resources+"/resources/nextclade_data/"+virus_dataset+"/reference.fasta --genemap="+path_to_resources+"/resources/nextclade_data/"+virus_dataset+"/genemap.gff --output-all="+output_dir+"/fasta_alignment_"+virus_dataset+" "+output_dir+"/temp_alignment_"+virus_dataset+"/*.fa" ,shell=True)

#one to do the agur
def create_phylogentic_tree(output_dir,runD):
    output_dir+="/"+runD
    #augur tree --method iqtree --alignment /home/ssh_user/output_nextali/nextalign.aligned.fasta 
    # --substitution-model GTR --nthreads 20 --output test_tree
    for virus_dataset in ["flu_h3n2_ha", "flu_h1n1pdm_ha"]:
        subprocess.run(". $CONDA_PREFIX/home/ssh_user/miniconda3/etc/profile.d/conda.sh && conda activate nextali && cd "+output_dir+"/fasta_alignment_"+virus_dataset+" && augur tree --method iqtree --substitution-model GTR --nthreads 20 --alignment "+output_dir+"/fasta_alignment_"+virus_dataset+"/nextalign.aligned.fasta --output "+output_dir+"/tree_"+virus_dataset+".nwk",shell=True)

        #copy alignment file out
        os.rename(output_dir+"/fasta_alignment_"+virus_dataset+"/nextalign.aligned.fasta",output_dir+"/"+runD+"_"+virus_dataset+"_aligned.fasta")

        #now clean up process
        subprocess.run("rm -r "+output_dir+"/temp_alignment_"+virus_dataset,shell=True)
        subprocess.run("rm -r "+output_dir+"/fasta_alignment_"+virus_dataset,shell=True)

def clean_run_files(pathTOirma,pathTOrundata,runD,pathToNextclade):

    subprocess.run("rm -r "+pathTOirma+runD,shell=True)
    subprocess.run("rm -r "+pathTOrundata+"/"+runD,shell=True)
    subprocess.run("rm -r "+pathToNextclade,shell=True)


if __name__ == "__main__":
    fastq_paths_dic={"2225102_060722_01" : "","2225196_060722_02" :"","2231833_060722_04" :" ","2225184_060722_03" :"",  "2229929_060722_05" :""}
    #move_fasta_files([*fastq_paths_dic],"/home/ssh_user/FLU_WGS_Sequencing/IRMA/","/home/ssh_user/FLU_WGS_Sequencing/results","060722")
    
    create_alignment_file("RESOURCES","/home/ssh_user/FLU_WGS_Sequencing/IRMA/",[*fastq_paths_dic],"/home/ssh_user/FLU_WGS_Sequencing/results","060722")

    create_phylogentic_tree("/home/ssh_user/FLU_WGS_Sequencing/results","060722")

   


