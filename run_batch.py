# python3 run_batch.py cd1 0 48 48

import sys, os

if len(sys.argv) < 1:
    print("Error")

# range of arg1 to arg2, eg run_batch.py 1 3 runs 1,2,3
# run_batch.py start_job end_job runtime_in_hours
for i in range(int(sys.argv[2]), int(sys.argv[3])+1):
    if os.path.exists("runqueue"):
        os.remove("runqueue")
    f = open("runqueue", "a")
    f.write("#!/bin/bash\n")
    f.write("#SBATCH -J run_queue_%s_%s\n" % (sys.argv[1], i))
    f.write("#SBATCH -p general\n")
    f.write("#SBATCH --mail-type=FAIL\n")
    f.write("#SBATCH --mail-user=sparpson@gmail.com\n")
    f.write("#SBATCH --mem=58G\n")
    f.write("#SBATCH --nodes=1\n")
    f.write("#SBATCH --ntasks=1\n")
    f.write("#SBATCH --cpus-per-task=48\n")
    f.write("#SBATCH --time=%s:00:00\n" % sys.argv[4])
    f.write("python3 /N/slate/jcandera/metaprod_projects/ml_project/manage.py runscript run_queue --script-args %s %s\n" % (sys.argv[1], i))
    f.close()
    os.system("sbatch runqueue")
