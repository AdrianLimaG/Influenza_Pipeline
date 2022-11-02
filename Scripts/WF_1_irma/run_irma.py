import time
from WF_1_irma.run_irma_helper import run_irma, move_results, sample_metrics

import json

#from merge_seq_data.merge_fastq_helper import find_samples

def irma_runner(dic_path,resource_path,results_dir):

    run_irma(dic_path, resource_path)
    time.sleep(10)

    move_results([*dic_path],resource_path,results_dir)

    metrics=sample_metrics([*dic_path],results_dir,resource_path)

    with open (resource_path+"/data/sample_metrics.json","w+") as j_dump :
        metrics = json.dumps(metrics)
        j_dump.write(metrics)



