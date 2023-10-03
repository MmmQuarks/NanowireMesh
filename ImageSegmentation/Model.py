import numpy as np
from shapely.geometry import Polygon, MultiPoint
import Algorithms as algos
_xy_to_idx, _idx_to_xy =  algos._xy_to_idx, algos._idx_to_xy

class Model(object):
	def __init__(self):
		super().__init__()
		self.img = None
		self.profiles = None
		self.polygon_points = []
		self.history = []
		self.track_history = True

	def cluster(self):
		# matplotlib assigns coordinates so that x is increasing from left to right 
		# and y is increasing from top to bottom. This means that 
		# x = column, # y = row
		new_img = np.zeros(self.img.shape)

		

	def highlight_nodes(self, *nodes):
		# making separate masks for each node to highlight
		try:
			bool_masks = [self.labeled_img.astype(np.double) == np.double(node) for node in nodes]

			# consecutive logical or on masks
			while len(bool_masks) > 1:
				bool_masks[0] = np.logical_or(bool_masks[0], bool_masks[1])
				del bool_masks[1]

			new_img =  np.where(
				bool_masks[0],
				2,
				np.where(
					self.labeled_img != 0,
					1,
					0
				)
			)
			return new_img
		except AttributeError:
			return None

	def clear_polygon_points(self):
		if self.polygon_points != []:
			print('Clearing polygon points. Previous values were \n{}'.format(
				'\n'.join(map(str,self.polygon_points))
			))
			self.polygon_points = []



	def segment(self, label):
		idx = np.argwhere(
			self.labeled_img == label
		)
		xy = [_idx_to_xy(tuple(el)) for el in idx]
		return MultiPoint(xy)

	def merge_segments(self, *segments):
		assert 0 not in segments, 'Cannot merge with a 0 labeled object'
		kept_label = min(segments)
		for label in segments:
			self.labeled_img[self.labeled_img == label] = kept_label
		return kept_label

	def new_segment_from_polygon_selection(self, *segments):
		polygon = Polygon(self.polygon_points)
		idx_to_relabel = []
		for label in segments:
			multipoint = polygon.intersection(
				self.segment(label)
			)
			for point in multipoint:
				xy = (point.x, point.y)
				idx_to_relabel.append(
					_xy_to_idx(xy)
				)
		new_label = self.labeled_img.max() + 1
		for idx in idx_to_relabel:
			row, col = idx
			row = int(row)
			col = int(col)
			self.labeled_img[row,col] = new_label
		return new_label