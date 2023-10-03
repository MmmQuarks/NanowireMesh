from pathlib import Path
import sys
import argparse

parser = argparse.ArgumentParser(
	prog = 'Max Power Dissipation Submit Generator',
	description = 'Modifies a template script into one or more customized scripts ready to be submitted with sbatch'
	)

parser.add_argument('target_files',
	nargs = '+',
	help = 'The target file(s) to be inserted into different version of the template script'
	)

if __name__ == '__main__':
	args = parser.parse_args(sys.argv[1:])
	for f in args.target_files:
		target_path = Path(f)
		job_name = 'max_power_dissipation_' + target_path.stem.replace('_nanowiremesh','')
		output = str(
				Path(
					'/home/amwt/TPV/2DNanowires/Data/2023-03-31AdamAgNW',
					job_name + '.out'
				)
			)
		error = str(
				Path(
					'/home/amwt/TPV/2DNanowires/Data/2023-03-31AdamAgNW',
					job_name + '.err'
				)
			)
		network_path = f

		template_path = '/home/amwt/TPV/2DNanowires/max_power_dissipation_template.py'

		template_text = Path(template_path).read_text()
		output_text = template_text.replace(
				'--job-name=template',
				'--job-name=' + job_name
			).replace(
				'--output=template.out',
				'--output=' + output
			).replace(
				'--error=template.err',
				'--error=' + error
			).replace(
				'<NETWORK_PATH>',
				network_path
			)
		output_path = Path(
				Path(template_path).parent,
				'submit_scripts',
				'submit_' + job_name + '.py'
				)
				
		template_path.replace(
			'template', job_name)
		with open(output_path, 'w') as output_file:
			output_file.write(output_text)
