from import_demo.import_demo_helper import demographics_import


def run_import_demo(runner_path,sample_hsn):

    import_demo = demographics_import(runner_path)

    import_demo.get_lims_demographics(sample_hsn)

    import_demo.format_lims_df()

    import_demo.merge_dfs()

    import_demo.format_dfs()

    import_demo.database_push()

    