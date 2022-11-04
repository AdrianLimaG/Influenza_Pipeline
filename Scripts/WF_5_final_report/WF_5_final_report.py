import os
import subprocess
from datetime import date



def create_final_report(run_date,path_to_nextclade,nextclade_hits,result_output):
    pass
    today = date.today().strftime("%b-%d-%Y")

    header=["Accession_Number","Patient_Name","Patient_DOB","Patient_Gender","Ordering_Facility", \
        "Specimen_Collection_Date","Race","Influenza_Type","Lineage_ID","Nextclade_Score","WGS_Run_Date","Report_Generation_Date"]
    #report_dir+"/"+date+"/"+date+"_demo.txt" this is the demo stuff
    #"Accession_Number","Patient_Name","Patient_DOB","Patient_Gender","Ordering_Facility","Specimen_Collection_Date","Race"

    #create report 
    report = open(result_output+"/"+run_date+"/"+run_date+"_report.tsv","w+")
    report.write("\t".join(header))

    #read in demo_file
    demo= open(result_output+"/"+run_date+"/"+run_date+"_demo.csv","r")

    #read_nextclade_file
    nextclade_file = open(path_to_nextclade+"/"+run_date+"/combined_nextclade.tsv")
    nextclade_dict={}
    for l in nextclade_file.readlines():
        l=l.strip()
        l=l.split("\t")
        nextclade_dict[l[0].split("_")[0]]=l[1],l[3]

    
    for line in demo.readlines() :
        #"Influenza_Type","Lineage_ID","Nextclade_Score","WGS_Run_Date","Report_Generation_Date"
        line = line.split(",")
        hsn = line[0]
                                #flu type                  #linageID            #score
        report_line = line + nextclade_hits[hsn][:] + [nextclade_dict[hsn][0]] + [nextclade_dict[hsn][1]]+ [run_date] + [today]


        report.write("\t".join(report_line))

    report.close()
    print("report created")

    





    #need to query the db for patient demographics or save demographical information in a earlier step


