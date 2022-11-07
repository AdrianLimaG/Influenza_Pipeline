from WF_3_nextclade.nextclade_helper import combine_nextclade_output, nextclade_data_obj, run_nextclade



def nextclade_runner(samples,results_dir,resources_path, nextclade_output_path):

    hits = run_nextclade(samples, results_dir, resources_path, nextclade_output_path)

    combined_file = combine_nextclade_output(nextclade_output_path,hits)
    
    nextclade_obj = nextclade_data_obj(resources_path)
    #print(dir(nextclade_obj))

    nextclade_obj.get_nextclade_dfs(combined_file)
    
    #commeting out to test and make sure everything works but this
    nextclade_obj.database_push()

    return hits


