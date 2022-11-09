import os
import subprocess
from datetime import date



def create_final_report(run_date,path_to_nextclade,nextclade_hits,result_output):
    pass
    today = date.today().strftime("%b-%d-%Y")
    hits={}

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

   
    for line in demo.readlines() :
        #"Influenza_Type","Lineage_ID","Nextclade_Score","WGS_Run_Date","Report_Generation_Date"
        line = line.strip().split(",")
        hsn = line[0]
                                #flu type                  #linageID            #score
        report_line = line + hits[hsn][:] + [nextclade_dict[hsn][0]] + [nextclade_dict[hsn][1]]+ [run_date] + [today]


        report.write("\t".join(report_line)+"\n")

    report.close()
    print("report created")


def move_fasta_files(fastq_samples,path_to_irma,result_output_dir,runD):

    os.mkdir(result_output_dir+"/"+runD+"/fasta_files")

    for sample in fastq_samples:
        temp_hsn = sample.split("_")[0]

        subprocess.run("cp "+path_to_irma+"/"+sample+"/amended_consensus/"+sample+"_combined.fasta "+result_output_dir+"/"+runD+"/fasta_files/"+temp_hsn+".fasta", shell=True)

    
    print("Fasta Files have been copied")



if __name__ == "__main__":
    fastq_paths_dic={"2225102_060722_01" : "","2225196_060722_02" :"","2231833_060722_04" :" ","2225184_060722_03" :"",  "2229929_060722_05" :""}
    
    move_fasta_files([*fastq_paths_dic],"/home/ssh_user/FLU_WGS_Sequencing/IRMA/","/home/ssh_user/FLU_WGS_Sequencing/results","060722")
    
   


