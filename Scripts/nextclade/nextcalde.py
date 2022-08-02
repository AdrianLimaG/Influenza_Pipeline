from nextclade.nextclade_helper import combine_nextclade_output, nextclade_data_obj, run_nextclade


#nextclade run  --in-order --input-dataset /home/ks_khel/data/flu_h3n2_ha/  --output-tsv /home/ks_khel/output/next.tsv --output-basename nextclade /home/ks_khel/Documents/GitHub/Infulenza_Pipeline/barcode01/amended_consensus/all_b1.fasta --output-all output


#needs to run nextcalde for each virus clade type if results is unsucessfull for all samples until 1 is succesfull or all have run
#starting with H3N1 H5N1 and the b yam ...
#need a function to run nextcalde
#need a function to handle nextclade outout

def nextclade_runner(samples,results_dir,resources_path, nextclade_output_path):

    hits = run_nextclade(samples, results_dir, resources_path, nextclade_output_path)
    #hits={'2225196': ['flu_h3n2_ha'], '2231833': ['flu_h3n2_ha'], '2225102': ['flu_h3n2_ha','flu_h1n1pdm_ha'], '2229929': ['flu_h3n2_ha'], '2225184': ['flu_h3n2_ha']}
        
    combined_file = combine_nextclade_output(nextclade_output_path,hits)

    nextclade_obj = nextclade_data_obj(resources_path)
    #print(dir(nextclade_obj))

    nextclade_obj.get_nextclade_dfs(combined_file)
    
    #commeting out to test and make sure everything works but this
    #nextclade_obj.database_push()

    return hits


