import sys
import time
from subprocess import Popen, PIPE, STDOUT


if __name__ == "__main__":
    print("="*50, "\nIN LINUX 'runner.py' SCRIPT\n" + "="*50)
    try:

        path_to_data = str(sys.argv[1])
        path_to_sampleSheet = str(sys.argv[2])
        # call to main script
        print("Path to Data "+path_to_data)
        print("Path to Sample Sheet "+path_to_sampleSheet)
        print("+"*50 + "\n$TNS_ADMIN BELOW\n" + "+"*50)
        p = Popen('export ORACLE_SID="KDHE_LIMS_CLOUD_PROD" && export ORACLE_HOME="/usr/lib/oracle/21/client64" && export TNS_ADMIN="/usr/lib/oracle/21/client64/network/admin/tnsnames.ora" && cat $TNS_ADMIN', executable="/usr/bin/bash", stdout=PIPE, stderr=STDOUT, shell=True)
        while True:
            line = p.stdout.readline()
            if not line:
                break
            print(line)
        
        p = Popen('xvfb-run export ORACLE_SID="KDHE_LIMS_CLOUD_PROD" && export ORACLE_HOME="/usr/lib/oracle/21/client64" && export TNS_ADMIN="/usr/lib/oracle/21/client64/network/admin/tnsnames.ora" && /home/ssh_user/miniconda3/bin/activate gen_workflow && /home/ssh_user/miniconda3/envs/gen_workflow/bin/python3.9 /home/ssh_user/Documents/GitHub/Infulenza_Pipeline/Scripts/flu_pipeline.py ' + path_to_data+' '+path_to_sampleSheet, executable="/usr/bin/sh", stdout=PIPE, stderr=STDOUT, shell=True)
        while True:
            line = p.stdout.readline()
            if not line:
                break
            print(line)
    except Exception as e:
        print(e)
        time.sleep(500)