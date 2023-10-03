import pandas as pd
import time
import numpy as np
import os

def main():
	df = pd.read_csv(
		'2023-03-31_AdamAgNW/NanowireDataIndex_scalebars.csv'
		)

	for idx, row in df.iterrows():
		if 'Image' in row.Type:
			if np.isnan(df.loc[idx, 'Scale Bar um']):
				os.system(
					'open {}'.format(
						row['New Path']
						)
					)
				ip = input('What is the scalebar size in pixels, then in um? Separate with a comma\n')
				pixels, um = [float(el) for el in ip.split(',')]
				df.loc[idx, 'Scale Bar um'] = um
				df.loc[idx, 'Scale Bar pix'] = pixels
				df.loc[idx, 'um per pix'] = um / pixels
				os.system(
					'osascript -e \"quit app \\\"Preview\\\"\"'
					)
		
	
				df.to_csv(
					'2023-03-31_AdamAgNW/NanowireDataIndex_scalebars.csv',
					index = False
					)
				time.sleep(0.5)
			else:
				print('already done')
				

if __name__ == '__main__':
	main()
