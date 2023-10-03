import numpy as np
import pandas as pd
from pathlib import Path
import ImageSegmentation.Algorithms as algos
import argparse
import sys
from angle_img_to_distance_matrix import main as distance_matrix
import pdb
from sklearn.cluster import DBSCAN as dbscan

def make_bins(df, nrows, ncols):
	"""
	Takes a dataframe corresponding do data from an image and gives and gives row/col bounds of subdivided image regions when image is divided into nrows rows and ncols columns
	args
		df = input dataframe with columns 'row' and 'col' at the very least
		nrows = the number of rows to cut it into
		ncols = the number of columns to cut it into
	returns
		row_bins, col_bins where these are the region boundaries for each sub-grid of the original image that the df has data from
	"""
		# subdividing the image
	minrow, mincol = df[['row','col']].min()
	maxrow, maxcol = df[['row','col']].max()
	row_bounds = np.linspace(
		minrow, 
		maxrow, 
		nrows + 1
	)
	col_bounds = np.linspace(
		mincol,
		maxcol,
		ncols + 1
	)
	row_bin_starts = row_bounds[:-1] * (1 - 0.125) # scaling these to make overlap happen
	col_bin_starts = col_bounds[:-1] * (1 - 0.125) 
	row_bin_starts = np.maximum(row_bin_starts, minrow).astype(int) # making sure we don't index out of bounds
	col_bin_starts = np.maximum(col_bin_starts, mincol).astype(int)
	row_bin_stops = row_bounds[1:] * 1.125
	col_bin_stops = col_bounds[1:] * 1.125
	row_bin_stops = np.minimum(row_bin_stops, maxrow).astype(int)
	col_bin_stops = np.minimum(col_bin_stops, maxcol).astype(int)

	row_bins = [el for el in zip(row_bin_starts, row_bin_stops)]
	col_bins = [el for el in zip(col_bin_starts, col_bin_stops)]
	return row_bins, col_bins


def main(**kwargs):
	"""
	Calculates a distance matrix based on a dataframe with columns x,y,angle,row,col
	where these correspond to data extracted from a radon_angle image
	args
		data 			str/Path to csv or pd.DataFrame object containing the data

		distance_threshold	only calculate full distance metric for points within this euclidean distance of each other

		save_patt		where to save the resulting distance matrix (if you want it saved at all). 
					This is taken as a prefix and appended with row/column information at save 
		
		nrows			Number of rows in the sub-grid of data

		ncols			Number of cols in the sub-grid of data
		
		eps			the maximum distance between two samples for one to be considered in the neighborhood of the other (i.e. for them to be direct neighbors) in dbscan

	returns
		None

	"""
	# allowing for data to either be a path or a pd.DataFrame file
	if type(kwargs['data']) in [str, type(Path())]:
		data_path = Path(kwargs['data'])
		df = pd.read_csv(data_path)
	elif type(kwargs['data']) == pd.DataFrame:
		df = kwargs['data']

	row_bins, col_bins = make_bins(df, kwargs['nrows'], kwargs['ncols'])

	for rb_idx in range(len(row_bins)):
		rb = row_bins[rb_idx]
		for cb_idx in range(len(col_bins)):
			print('Image Region {}/{}'.format(
				rb_idx * len(col_bins) + cb_idx,
				len(row_bins) * len(col_bins) - 1
				)
			)
			cb = col_bins[cb_idx]
			rcdf = df[
				(df.row >= rb[0]) & (df.row <= rb[1]) &
				(df.col >= cb[0]) & (df.col <= cb[1])
				].reset_index(drop = False)
			file_pattern = '{}_dist_mat_row_{}_col_{}'.format(
				kwargs['save_patt'],
				int(rb_idx),
				int(cb_idx)
				)

			dist_mat = distance_matrix(
				data = rcdf,
				distance_threshold = kwargs['distance_threshold']
				)
			print('Clustering')
			model = dbscan(
				metric = 'precomputed',
				eps = kwargs['eps']
				)
			rcdf['labels'] = model.fit_predict(dist_mat)
			rcdf.to_csv(
				file_pattern + '.csv'
				)
			



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
		'--save-patt',
		help = 'the pattern to follow (including a path) when naming chunked dbscan data'
		)
	parser.add_argument(
		'--nrows',
		type = int,
		help = 'The number of rows in the sub-grid of data',
		default = 1
		)
	parser.add_argument(
		'--ncols',
		type = int,
		help = 'The number of cols of the sub-grid of data',
		default = 1
		)
	parser.add_argument(
		'--eps',
		type = float,
		help = 'the maximum distance between two samples for one to be considered in the neighborhood of the other (i.e. for them to be direct neighbors) in dbscan',
		default = 5
		)

	args = parser.parse_args(sys.argv[1:])
	main(
		**vars(args)
	)

