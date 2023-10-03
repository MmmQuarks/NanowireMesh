import tkinter as tk
import pdb
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.collections import LineCollection
import matplotlib as mpl
from numpy import array


class ImageView(tk.Frame):
	def __init__(self, controller = None,**kwargs):
		assert controller is not None, 'Must assign a controller when instantiating an ImageView'
		super().__init__(**kwargs)
		self.controller = controller
		self._canvas = None
		self._ax = None
		self.mpl_img = None
		self.mpl_plot = None
		self.mpl_polygon = None
		self.first_internal_plot = True
	
	@property
	def img(self):
		return self.controller.model.img

	@property
	def ax(self):
		if self._ax is None:
			fig = Figure()
			fig.set_facecolor('black')
			self._ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
			self._ax.set_axis_off()
			self._ax.set_aspect('equal')
			self._ax.set_facecolor('black')
		return self._ax

	@property
	def canvas(self):
		if self._canvas is None:
			# fig = Figure()
			# self.ax = fig.add_axes([0.05, 0.05, .9, .9])
			# self.ax.set_axis_off()
			self._canvas = FigureCanvasTkAgg(
				figure = self.ax.get_figure(),
				master = self
			)
			cid = self._canvas.mpl_connect(
				'button_press_event',
				self.dbl_click_function
			)
			toolbar = NavigationToolbar2Tk(
				self.canvas,
				self,
				pack_toolbar = False
			)
			toolbar.pack(side = tk.BOTTOM)
			self.canvas.get_tk_widget().pack(
				side = tk.TOP,
				fill = tk.BOTH,
				expand = True
			)
		return self._canvas
	
	def set_img(self):
		'''
		to change the image being seen, set the img attribute of controller 
		'''
		if self.mpl_img == None:
			# fig = Figure()
			# self.ax = fig.add_axes([0.05, 0.05, .9, .9])
			# self.ax.set_axis_off()

			# confirm canvas exists
			self.canvas
			# make img
			self.mpl_img = self.ax.imshow(self.controller.model.img, cmap = 'hot')
			# self.mpl_plot = None
			# self.canvas = FigureCanvasTkAgg(
				# figure = fig,
				# master = self
			# )
			# cid = self.canvas.mpl_connect(
			# 	'button_press_event',
			# 	self.dbl_click_function
			# )
			self.canvas.draw()
			# toolbar = NavigationToolbar2Tk(
			# 	self.canvas,
			# 	self,
			# 	pack_toolbar = False
			# )
			# toolbar.pack(side = tk.BOTTOM)
			# self.canvas.get_tk_widget().pack(
			# 	side = tk.TOP,
			# 	fill = tk.BOTH,
			# 	expand = True
			# )

		else:
			self.mpl_img.set_data(self.controller.model.img)
			self.mpl_img.set_clim(
				vmin = self.controller.model.img.min(),
				vmax = self.controller.model.img.max()
			)
			self.canvas.draw()

	def set_plot(self):
		g = self.controller.model.g
		nodes_selected = self.controller.nodes_selected

		#make sure canvas and axis objects exist
		self.canvas
		if nodes_selected is None:
			sg = g
		else:
			sg = g.subgraph(nodes_selected)
		

		segment_coords = [list(sg.nodes[node]['shape'].coords[:]) for node in sg]
		try: # if mpl_segments exists, modify it
			self.mpl_segments.set_segments(
				segment_coords
			)
		except AttributeError: # if mpl_segments not exist yet, create it 
			segment_collection = LineCollection(
				segment_coords,
				linewidths = 1,
				colors = 'blue'
			)
			self.mpl_segments = self.ax.add_collection(segment_collection)

		# now do similar thing for the resistors (but only if they are there)
		markerSize = 1.5 * mpl.rcParams['lines.markersize']**2

		# all contact edges shown here
		contacts = [edge for edge in g.edges if g.edges[edge]['resistanceType'] == 'cont']
		if contacts:
			# only contacts in the subgraph of selected nodes will have nonzero size
			sizes = [markerSize if edge in sg.edges else 0 for edge in contacts]
			try:
				# try to set sizes (will throw errorw if plot not exist yet)
				self.mpl_contact_resistors.set_sizes(sizes)#.remove()
			except AttributeError:
				# if plot not exist yet, make plot
				contact_resistor_coords = array(
					[g.edges[edge]['shape'].coords[0] for edge in contacts]
				)
				self.mpl_contact_resistors = self.ax.scatter(
					x = contact_resistor_coords[:,0],
					y = contact_resistor_coords[:,1],
					marker = 'o',
					c = 'cyan',
					s = sizes
				)

		# all internal edges shown here
		internals = [edge for edge in g.edges if g.edges[edge]['resistanceType'] == 'int']
		if internals:
			# only internal edges in the subgraph of selected nodes will have nonzero size
			sizes = [markerSize if edge in sg.edges else 0 for edge in internals]
			try:
				self.mpl_internal_resistors.set_sizes(sizes)#.remove()
				# try to set sizes (will throw error if plot not exist yet)
			except AttributeError:
				# if plot not exist, make the plot
				internal_resistor_coords = array(
					[g.edges[edge]['shape'].coords[0] for edge in internals]
				)
				self.mpl_internal_resistors = self.ax.scatter(
					x = internal_resistor_coords[:,0],
					y = internal_resistor_coords[:,1],
					marker = '*',
					c = 'lime',
					s = sizes
				)


		
		# redraw the plot
		self.ax.set_facecolor('black')
		self.canvas.draw()




	def set_controller(self, controller):
		self.controller = controller

	def dbl_click_function(self, event):
		if event.dblclick:
			# processing the click data
			displayToDataCoordTransformer = self.ax.transData.inverted()
			x,y = displayToDataCoordTransformer.transform((event.x, event.y))
			if self.controller.click_mode.get() == 'Select Polygon':
				self.controller.model.polygon_points.append(
					(x, y)
				)
				self.controller.show_polygon()
			elif self.controller.click_mode.get() == 'Click to Highlight':
				self.controller.add_highlight_to_nearest(x,y)
			else:
				print('xdata = {}, ydata = {}'.format(event.xdata, event.ydata))
				displayToData = self.ax.transData.inverted()
				print('converted from display coordinates')
				print(displayToData.transform((event.x, event.y)))