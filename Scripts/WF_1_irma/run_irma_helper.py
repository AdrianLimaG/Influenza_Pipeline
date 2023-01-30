import os
import subprocess
import time



def run_irma(dic_path_data,path_to_irma):
    #this Function runs IRMA 
    #inpput should be a string that is that path of current repository to which i will add proper extion to find irm
    #2 input should be an array of string which lead to path of combined fastq.gz
    #in future should probaly be a dictonary with HSN followed by path
    
    irma_path = path_to_irma+"/resources/flu-amd/"
    
    irma_command= irma_path+"IRMA FLU-pgm "
    #output path should be a variable passed down based on parameter of function
    for key in dic_path_data:
        irma_command= irma_path+"IRMA FLU-pgm "
        irma_command+=dic_path_data[key]+" "+key
        #print(irma_command)
        subprocess.run(irma_command,shell=True)
        time.sleep(10)


def move_results(samples,irma_output_dir,results_dir):
    #IRMA only outputs in IRMA's working directory
    #This function moves the results to there permant home
    #checking is date folder has been created if not making it
    if not os.path.exists(results_dir):
        subprocess.run("mkdir "+results_dir,shell=True)

    for sample in samples:
        subprocess.run("mv "+irma_output_dir+"/"+sample+" "+results_dir+"/"+sample, shell=True)

def sample_metrics(list_samples,sample_path,resource_path):
    #This functions calculates WGS depth + coverage an anverage across all aligned proteins
    #And Captures this for each protein and writes it to a json file
    samtools_path = resource_path+"/resources/samtools/bin/samtools coverage "
   # coverage={} #return vaule containing key being sample_name
    c2=[]
    #coverage = { "Sample1": {"WGS_Coverage": xx, "WGS_Depth": xx, "HA_Coverage": xx , "HA_Depth" : xx}     }
    
    for sample in list_samples:
        
        temp_depth=[]
        temp_coverage=[]
        hsn= sample.split("_")[0]
        temp_dict ={"hsn": hsn}
        #coverage[hsn] = {}
        #[]
        for protein in ["A_HA_H3.bam","A_MP.bam","A_NA_N2.bam","A_NP.bam","A_NS.bam","A_PA.bam","A_PB1.bam","A_PB2.bam"] :

            if os.path.exists(sample_path+"/"+sample+"/"+protein):

                out_put=subprocess.run(samtools_path+sample_path+"/"+sample+"/"+protein, capture_output=True, text=True, shell=True)
                print(out_put)
                data= out_put.stdout.split("\n")[1].split("\t")
                temp_dict[protein[2:-4]+"_avg_depth"] = round(float(data[6]),2)
                temp_dict[protein[2:-4]+"_coverage"] = round(float(data[5]),2)
                #coverage[hsn][protein[2:-4]+"_avg_depth"] = int(round(float(data[6])))
                #coverage[hsn][protein[2:-4]+"_coverage"] = int(data[5])
                temp_depth.append(round(float(data[6]),2))
                temp_coverage.append(round(float(data[5]),2))

            else:
                print("sample "+sample+" has no bam file for protein "+ protein)  #should be written to a file  
                temp_dict[protein[2:-4]+"_avg_depth"] = 0
                temp_dict[protein[2:-4]+"_coverage"] = 0
               # coverage[hsn][protein[2:-4]+"_avg_depth"] = 0
               # coverage[hsn][protein[2:-4]+"_coverage"] = 0
                temp_depth.append(0)
                temp_coverage.append(0)
                
        temp_dict["percent_cvg"]= sum(temp_coverage)/len(temp_coverage)
        temp_dict["avg_depth"]= sum(temp_depth)/len(temp_depth)
        #coverage[hsn]["percent_cvg"]= sum(temp_coverage)/len(temp_coverage)
        #coverage[hsn]["avg_depth"]= sum(temp_depth)/len(temp_depth)
        c2.append(temp_dict)
        
    
    #print(coverage)
    return c2

            


    

'''if __name__ == "__main__":
    d= {'2225102': '/home/ks_khel/Desktop/FLU_DATA/barcode01/2225102_combined.fastq.gz', '2229929': '/home/ks_khel/Desktop/FLU_DATA/barcode05/2229929_combined.fastq.gz', '2225184': '/home/ks_khel/Desktop/FLU_DATA/barcode03/2225184_combined.fastq.gz'}
    run_irma(d,"/home/ks_khel/Documents/GitHub/Infulenza_Pipeline")
    move_results([*d],"/home/ks_khel/Documents/GitHub/Infulenza_Pipeline","/home/ks_khel/Desktop/RES/")'''
''