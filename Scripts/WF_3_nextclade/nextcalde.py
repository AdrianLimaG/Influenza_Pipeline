from WF_3_nextclade.nextclade_helper import combine_nextclade_output, nextclade_data_obj, run_nextclade



def nextclade_runner(samples,results_dir,resources_path, nextclade_output_path):

    hits = run_nextclade(samples, results_dir, resources_path, nextclade_output_path)

    #hits={'2225196': ['flu_h3n2_ha'], '2231833': ['flu_h3n2_ha'], '2225102': ['flu_h3n2_ha','flu_h1n1pdm_ha'], '2229929': ['flu_h3n2_ha'], '2225184': ['flu_h3n2_ha']}
    #print(hits)    

    combined_file = combine_nextclade_output(nextclade_output_path,hits)
    #combined_file = "/home/ssh_user/FLU_WGS_Sequencing/Nextclade/combined_nextclade.tsv"
    
    nextclade_obj = nextclade_data_obj(resources_path)
    #print(dir(nextclade_obj))

    nextclade_obj.get_nextclade_dfs(combined_file)
    
    #commeting out to test and make sure everything works but this
    nextclade_obj.database_push()

    return hits


