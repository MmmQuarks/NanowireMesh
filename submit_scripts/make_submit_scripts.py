#python script to write submit scripts for randomized networks analysis

preamble = "#!/bin/bash\n"
preamble = preamble + "#\n"
preamble = preamble + "#SBATCH --job-name=spice_analyze_network\n"
preamble = preamble + "#SBATCH --ntasks=1\n" + "#SBATCH --mem-per-cpu=4000\n" + "#SBATCH -t 30:00:00\n\n"

for taskId in range(1,21):
	for iterId in range(1,6):
		pyCommand = "srun python ~/TPV/3D_nanowires/pythonSpiceManager.py " + str(taskId) + " " + str(iterId)
		fileName = "submit_" + str(taskId) + "_" + str(iterId)
		f = open(fileName, "w")
		f.write(preamble + pyCommand)
		f.close()
