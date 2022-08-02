import time
from irma.run_irma_helper import run_irma
from irma.run_irma_helper import move_results

#from merge_seq_data.merge_fastq_helper import find_samples

def irma_runner(dic_path,irma_path,results_dir):

    run_irma(dic_path, irma_path)
    time.sleep(10)
    move_results([*dic_path],irma_path,results_dir)

