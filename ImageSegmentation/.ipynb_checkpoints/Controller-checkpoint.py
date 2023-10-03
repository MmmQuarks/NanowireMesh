import tkinter as tk
import pdb
import pickle
import Algorithms as algos
from to_nanowire_mesh import to_nanowire_mesh
import numpy as np
from matplotlib.lines import Line2D
from skimage import io
from random import sample
import os


class Controller(object):
	def __init__(self, model = None):
		super().__init__()
		self.model = model
		self.polygon_mode = tk.IntVar()
		self.add_highlight_to_nearest_mode = tk.IntVar()
		self.click_mode = tk.StringVar()
	
	@property
	def img(self):
		return self.model.img

	@img.setter
	def img(self, new_img):
		if new_img is not None:
			self.model.img = new_img
			self.image_view.set_img()
		else:
			print('Image is None. No action taken.')

	@property
	def nodes_selected(self):
		if self.action_view.ent_nodes_selected.get() == '':
			return None
		selected_string = str.strip(
			self.action_view.ent_nodes_selected.get()
		)
		selected = selected_string.split(',')
		selected = [float(el) for el in selected]
		return selected

	def log(self, text):
		self.action_view.txt_log.insert(tk.INSERT, text + 2 * '\n')

	def set_nodes_selected(self, *nodes):
		string_nodes = map(str, nodes)
		selected_string = ','.join(string_nodes)
		self.action_view.ent_nodes_selected.delete(0, tk.END)
		self.action_view.ent_nodes_selected.insert(0, selected_string)


	def set_image_view(self,image_view):
		self.image_view = image_view

	def set_action_view(self, action_view):
		self.action_view = action_view

	def load_profiles_and_img(self, file_path = None):
		if file_path is None:
			file_path = tk.filedialog.askopenfilename(
				initialdir = os.getcwd(),
				title = 'Select Pickle'
			)
		with open(file_path, 'rb') as profiles_file:
			pickled_data = pickle.load(profiles_file)
			self.model.profiles = pickled_data['profiles']
			self.img = pickled_data['img']
			self.model.original_img = pickled_data['img']

		labels = algos.cluster_profiles(self.model.profiles)
		pixel_idx = list(self.model.profiles.keys())
		assert len(labels) == len(pixel_idx)
		new_img = np.zeros(self.img.shape)
		for n in range(len(labels)):
			row,col = pixel_idx[n]
			new_img[row,col] = labels[n]
		self.model.labeled_img = new_img
		self.img = new_img

	def load_graph(self, file_path = None):
		if file_path is None:
			file_path = tk.filedialog.askopenfilename(
				initialdir = os.getcwd(),
				title = 'Select Pickle'
			)
		with open(file_path, 'rb') as graph_file:
			self.model.g = pickle.load(graph_file)
		self.set_plot()
		self.print_parallel_edges()

	def save_graph(self, file_path = None):
		if file_path is None:
			file_path = tk.filedialog.asksaveasfilename(
				initialdir = os.getcwd(),
				defaultextension = '.pickle'
			)
		with open(file_path, 'wb') as graph_file:
			pickle.dump(self.model.g, graph_file)


	def save_img(self):
		file_name = tk.filedialog.asksaveasfilename(
			initialfile = 'Untitled.png',
			initialdir = os.getcwd(),
			defaultextension = '.png'
		)
		img = algos.gray_to_rgb(self.img)
		io.imsave(
			file_name,
			img
		)
		self.model.saveimg = img

	def load_labeled_img(self, file_path = None):
		if file_path is None:
			file_path = tk.filedialog.askopenfilename(
				initialdir = '/',
				title = 'Select image for import'
			)
		labeled_img = algos.rgb_to_gray(
			io.imread(file_path)
		)

		self.img = self.model.original_img = self.model.labeled_img = labeled_img


	def home(self):
		self.img = self.model.labeled_img
		self.action_view.ent_nodes_selected.delete(0, tk.END)
		self.set_plot()

	def single_node_view(self):
		nodes = np.unique(self.model.labeled_img)
		self.set_nodes_selected(nodes[0])
		self.img = self.model.highlight_nodes(
			self.nodes_selected[0]
		)


	def go_to_node(self):
		self.img = self.model.highlight_nodes(
			*self.nodes_selected
		)
		self.image_view.set_plot()

	def next_node(self):
		# if we have multiple nodes that are selected currently, we go to the second node currently selected
		if len(self.nodes_selected) > 1:
			self.set_nodes_selected(
				self.nodes_selected[1]
			)
		else: # if we only have one node selected, we go to the next one in the sorted list of nodes
			nodes = sorted(list(np.unique(self.model.labeled_img)))
			current_idx = nodes.index(self.nodes_selected[0])
			next_idx = (current_idx + 1) % len(nodes)
			self.set_nodes_selected(
				nodes[next_idx]
			)
		
		# changing the img
		self.go_to_node()


	def prev_node(self):
		# if we have multiple nodes that are selected currently, we go to the second node currently selected
		if len(self.nodes_selected) > 1:
			self.set_nodes_selected(
				self.nodes_selected[0]
			)
		else: # if we only have one node selected, we go to the next one in the sorted list of nodes
			nodes = sorted(list(np.unique(self.model.labeled_img)))
			current_idx = nodes.index(self.nodes_selected[0])
			prev_idx = current_idx - 1
			self.set_nodes_selected(
				nodes[prev_idx]
			)

		# changing the img
		self.go_to_node()

	def add_highlight_to_nearest(self, x, y):
		labeled_img = self.model.labeled_img
		nearest_idx = algos.nearest_positive_idx(labeled_img, x, y)
		label = labeled_img[nearest_idx[0], nearest_idx[1]]
		selected = self.nodes_selected
		if label not in selected:
			selected.append(label)
			self.set_nodes_selected(*selected)
			self.go_to_node()


	def show_polygon(self):
		self.unshow_polygon()
		polygon_points = np.array(self.model.polygon_points)
		self.image_view.line = Line2D(
			polygon_points[:,0],
			polygon_points[:,1]
		)
		self.image_view.mpl_line = self.image_view.ax.add_line(
			self.image_view.line
			)
		self.image_view.canvas.draw()

	def unshow_polygon(self):
		try:
			self.image_view.mpl_line.remove()
			del self.image_view.mpl_line
			self.image_view.canvas.draw()
		except AttributeError:
			#nothing to delete
			pass

	def delete_polygon(self):
		self.unshow_polygon()
		self.model.polygon_points = []

	def delete_last_polygon_point(self):
		self.unshow_polygon()
		del self.model.polygon_points[-1]
		self.show_polygon()


	def merge_segments(self):
		if len(self.nodes_selected) < 2:
			pass
		else:
			kept_label = self.model.merge_segments(
				*self.nodes_selected
			)
			self.set_nodes_selected(
				kept_label
			)	
			self.go_to_node()

	def new_segment_from_polygon_selection(self):
		new_label = self.model.new_segment_from_polygon_selection(
			*self.nodes_selected
		)
		self.set_nodes_selected(
			new_label
		)
		self.go_to_node()
		self.delete_polygon()

	# def check_intersections(self):
	# 	g = self.model.g
	# 	# find all pairs of nodes that have more than one edge between them
	# 	multiple_edges = set(
	# 		[(e[0], e[1]) for e in g.edges if g.number_of_edges(e[0], e[1]) > 1]
	# 	)

	# 	# now for each pair of these plot it and highlight the nodes in question



	def set_plot(self):
		self.image_view.set_plot()

	def labeled_img_to_graph(self):
		self.model.g = algos.labeled_img_to_graph(
			img = self.model.labeled_img,
			add_electrodes = True
		)

	def to_nanowire_mesh(self):
		nwDiam = float(
			self.action_view.ent_nwDiam.get()
		)
		initialTemp = float(
			self.action_view.ent_initialTemp.get()
		)
		pixels = float(
			self.action_view.ent_scale_bar_pixels.get()
		)
		um = float(
			self.action_view.ent_scale_bar_um.get()
		)
		self.log('Making Nanowire Mesh, including internal resistances, with')
		self.log('nwDiam = {}'.format(nwDiam))
		self.log('initialTemp = {}'.format(initialTemp))
		self.log('pixels per um = {}'.format(pixels / um))
		self.log('removeNonPercolationg = {}'.format(False))
		self.log('contact resistance distribution = bellew')


		self.model.g = to_nanowire_mesh(
			self.model.g,
			pixels = pixels,
			um = um,
			nwDiam = nwDiam,
			initialTemp = initialTemp
			)
		self.set_plot()

	def print_parallel_edges(self):
		pe = algos.parallel_edges(self.model.g)
		msg = '\n'.join(
			['Nodes {} have edges {}'.format(key,val) for key,val in pe.items()]
		)
		msg = msg
		self.log(msg)