import numpy as np
import pandas as pd
from pathlib import Path
import ImageSegmentation.Algorithms as algos
import argparse
import sys



def main(**kwargs):
	"""
	Calculates a distance matrix based on a dataframe with columns x,y,angle,row,col
	where these correspond to data extracted from a radon_angle image
	args
		data 			str/Path to csv or pd.DataFrame object containing the data

		distance_threshold	only calculate full distance metric for points within this euclidean distance of each other

		save_dist_mat		path for saving the output distance matrix
		
	returns
		if save_dist_mat is not specified,
			returns dist_mat
		else:
			saves dist_mat and returns none

	"""
	if type(kwargs['data']) in [str, type(Path())]:
		data_path = Path(kwargs['data'])
		df = pd.read_csv(data_path)
	elif type(kwargs['data']) == pd.DataFrame:
		df = kwargs['data']
	distance_threshold = kwargs['distance_threshold']
	dist_mat = np.ones(
		(len(df), len(df))
		) * 1E5 # so that the distance will be very large between points
		# separated by more than the euclidean distance threshold
	for idx, row in df.iterrows():
		if idx % 1E3 == 0 or idx == len(df) - 1:
			print('Distance Matrix chunk {}/{}'.format(idx, len(df) - 1))
	
		# we consider p1 to be the left point and p2 to be the right point
		# hence the name rdf (which stands for right df which contains all the p2)
		rdf = df[df.index > idx].copy()
		p1 = row[['x','y']].values
		p2 = rdf[['x','y']].values
		rdf['distance'] = np.linalg.norm(
			p1 - p2,
			axis = 1
			)
		# n stands for neighbor df
		ndf = rdf[rdf.distance <= distance_threshold].copy()
		if len(ndf) > 0:
			ndf['angle_diff'] = algos.line_angle_difference(
				row.angle,
				ndf.angle
				)
			dist_mat_row_idx = [idx] * len(ndf)
			dist_mat_col_idx = ndf.index.values
			dist_mat[dist_mat_row_idx, dist_mat_col_idx] = ndf.distance.values + ndf.angle_diff.values
	

	# note that, if we are executing this as a script, the first block here will always be executed because it is required by argparse
	if 'save_dist_mat' in kwargs.keys():
		print('Attempting to save the distance matrix')
		np.save(
			kwargs['save_dist_mat'],
			dist_mat
			)
	return dist_mat

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'--data',
		help = 'Path to csv containing columns x,y,angle,row,col'
		)
	parser.add_argument(
		'--distance-threshold',
		type = float,
		default = 5,
		help = 'The distance between pixels, in distance units of the dataset, for them to be possibly considered direct neighbors'
		)
	parser.add_argument(
		'--save-dist-mat',
		help = 'where to save the distance matrix',
		required = True
		)
	args = parser.parse_args(sys.argv[1:])

	main(**vars(args))


