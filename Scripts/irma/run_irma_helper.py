import os
import subprocess
import time



def run_irma(dic_path_data,path_to_irma):
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

    for sample in samples:
        subprocess.run("mv "+irma_output_dir+"/"+sample+" "+results_dir, shell=True)

'''if __name__ == "__main__":
    d= {'2225102': '/home/ks_khel/Desktop/FLU_DATA/barcode01/2225102_combined.fastq.gz', '2229929': '/home/ks_khel/Desktop/FLU_DATA/barcode05/2229929_combined.fastq.gz', '2225184': '/home/ks_khel/Desktop/FLU_DATA/barcode03/2225184_combined.fastq.gz'}
    run_irma(d,"/home/ks_khel/Documents/GitHub/Infulenza_Pipeline")
    move_results([*d],"/home/ks_khel/Documents/GitHub/Infulenza_Pipeline","/home/ks_khel/Desktop/RES/")'''
