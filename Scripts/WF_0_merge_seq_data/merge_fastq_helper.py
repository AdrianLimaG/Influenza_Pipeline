
import os
import subprocess


def find_samples(path_to_data,sample_sheet_p):
    folders = [x[0] for x in os.walk(path_to_data)]

    combine_file_paths={}
    hsns = read_sample_sheet(sample_sheet_p)
    #skipping parent directory
    for sample_folder_path in folders[1:]:
        combine_file_paths.update(concat_fastq(sample_folder_path,hsns))

    return combine_file_paths

def concat_fastq(path_to_fastq,hsn_dict):
    fastq = os.listdir(path_to_fastq)
    fastq_s = (" "+path_to_fastq+"/").join(['',*fastq])

    barcode = path_to_fastq.split("/")[-1]
    
    #need to map barcode to HSN now as this will make it very easy
    #findposiion
    hsn= hsn_dict[barcode]

    outputfilepath= path_to_fastq+"/"+hsn+"_combined.fastq.gz"
    #print("cat "+fastq_s+" > "+outputfilepath)
    output= subprocess.run("cat "+fastq_s+" > "+outputfilepath,shell=True)

    return {hsn:outputfilepath}


def read_sample_sheet(sample_sheet_path):

    sample_sheet = open(sample_sheet_path,"r")
    sample_sheet_date = sample_sheet_path.split("/")[-1][:-4]  #geting the date from the samplesheet file name
    lines=sample_sheet.readlines()
    hsn={}

    for line in lines[1:]:
        l=line.strip().split(",")
        #l[0] HSN
        #l[1] plate pos
        #l[2] barcode
        if int(l[2]) < 10 :
            l[2]="0"+l[2]
        hsn["barcode"+l[2]]=l[0]+"_"+sample_sheet_date+"_"+l[2]
    sample_sheet.close()

    return hsn



if __name__ == "__main__":
    
    sample_sheet_p ="/home/ks_khel/Desktop/FLU_DATA/060722.csv" #will need a way to pass this differently
    read_sample_sheet(sample_sheet_p)

