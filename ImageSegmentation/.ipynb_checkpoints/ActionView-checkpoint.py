import tkinter as tk
class ActionView(tk.Frame):
	def __init__(self, controller = None, **kwargs):
		assert controller is not None, 'Must assign a controller when instantiating an ActionView'
		super().__init__(**kwargs)
		self.controller = controller

		self.columnconfigure([0,2], weight = 1, pad = 5, minsize = 25)
		self.columnconfigure(1, weight = 2, pad = 5, minsize = 100)
		self.rowconfigure(2, minsize = 50)

		# Home Button
		btn_home = tk.Button(
			master = self,
			text = 'Home',
			command = self.controller.home
		)
		btn_home.grid(row = 2, column = 0, sticky = 'nsew')


		# save graph
		btn_save_graph = tk.Button(
			master = self,
			text = "Save Graph",
			command = self.controller.save_graph
		)
		btn_save_graph.grid(row = 2, column = 1)

		# Load profiles and image button
		btn_load_profiles_and_img = tk.Button(
			master = self,
			text = 'Load Profiles\nand Image',
			command = self.controller.load_profiles_and_img
		)
		btn_load_profiles_and_img.grid(row = 1, column = 1)

		# load graph pickle
		btn_load_graph = tk.Button(
			master = self,
			text = 'Load Graph',
			command = self.controller.load_graph
		)
		btn_load_graph.grid(row = 1, column = 0)

		# Load labeled image button
		btn_load_labeled_img = tk.Button(
			master = self,
			text = 'Load Labeled Image',
			command = self.controller.load_labeled_img
		)
		btn_load_labeled_img.grid(row = 1, column = 2)

		# save image button
		btn_save_img = tk.Button(
			master = self,
			text = 'Save Image',
			command = self.controller.save_img
		)
		btn_save_img.grid(row = 2, column = 2, sticky = 'nsew')

		# single node view button
		btn_single_node_view = tk.Button(
			master = self,
			text = 'Single Node View',
			command = self.controller.single_node_view
		)
		btn_single_node_view.grid(row = 3, column = 1)

		# nodes selected text box
		self.ent_nodes_selected = tk.Entry(master = self)
		self.ent_nodes_selected.grid(row = 4, column = 1)
		self.ent_nodes_selected.bind('<Return>', lambda event : self.controller.go_to_node())

		next_prev_node_frame = tk.Frame(master = self)
		next_prev_node_frame.grid(row = 4, column = 2)

		# show next/previous node buttons
		btn_next_node = tk.Button(
			master = next_prev_node_frame,
			text = '-->',
			command = self.controller.next_node
		)
		btn_prev_node = tk.Button(
			master = next_prev_node_frame,
			text = '<--',
			command = self.controller.prev_node
		)
		btn_prev_node.grid(row = 0, column = 0)
		btn_next_node.grid(row = 0, column = 1)

		# go to node(s) displayed in ent_nodes_selected
		btn_go_to_node = tk.Button(
			master = self,
			text = 'Go',
			command = self.controller.go_to_node
		)
		btn_go_to_node.grid(row = 4, column = 0)
		self.rowconfigure(5, minsize = 50)

		# click mode radio buttons
		radio_none_mode = tk.Radiobutton(
			master = self,
			text = 'Inactive',
			variable = self.controller.click_mode,
			value = None,
		)
		radio_none_mode.grid(row = 6, column = 1)
		radio_polygon_mode = tk.Radiobutton(
			master = self,
			text = 'Select Polygon',
			variable = self.controller.click_mode,
			value = 'Select Polygon',
			command = self.controller.model.clear_polygon_points
		)
		radio_polygon_mode.grid(row = 7, column = 1)

		radio_add_highlight_to_nearest = tk.Radiobutton(
			master = self,
			text = 'Click to Highlight',
			variable = self.controller.click_mode,
			value = 'Click to Highlight'
		)
		radio_add_highlight_to_nearest.grid(row = 8, column = 1)


		# hide polygon button
		btn_unshow_polygon = tk.Button(
			master = self,
			text = 'Unshow Polygon',
			command = self.controller.unshow_polygon
		)
		btn_unshow_polygon.grid(row = 6, column = 0)

		# delete last polygon point button
		btn_delete_last_polygon_point = tk.Button(
			master = self,
			text = 'Delete Last\nPolygon Point',
			command = self.controller.delete_last_polygon_point
		)
		btn_delete_polygon = tk.Button(
			master = self,
			text = 'Delete Polygon',
			command = self.controller.delete_polygon
		)

		btn_delete_last_polygon_point.grid(row = 6, column = 2)
		btn_delete_polygon.grid(row = 7, column = 2)

		# new segment from polygon selection button
		btn_new_segment_from_polygon_selection = tk.Button(
			master = self,
			text = 'New Segment\nfrom Polygon',
			command = self.controller.new_segment_from_polygon_selection
		)
		# merge segments in text box button
		btn_merge_segments = tk.Button(
			master = self,
			text = 'Merge Segments\nin Text Box',
			command = self.controller.merge_segments
		)
		btn_new_segment_from_polygon_selection.grid(row = 8, column = 0)
		btn_merge_segments.grid(row = 8, column = 2)


		lbl_network_actions = tk.Label(
			master = self,
			text = '\nNetwork Actions\n',
			font = ('Arial' , 25)
		)
		lbl_network_actions.grid(row = 9, column = 1)

		btn_find_intersections = tk.Button(
			master = self,
			text = 'Labeled Image to Graph \n(no internal resistance)',
			command = self.controller.labeled_img_to_graph
		)
		lbl_initialTemp = tk.Label(
			master = self,
			text = 'initialTemp [K]'
		)
		ent_initialTemp = tk.Entry(
			master = self
			)
		lbl_nwDiam = tk.Label(
			master = self,
			text = 'nwDiam [um]'
		)
		ent_nwDiam = tk.Entry(
			master = self
		)

		lbl_initialTemp.grid(row = 10, column = 0)
		ent_initialTemp.grid(row = 11, column = 0)
		lbl_nwDiam.grid(row = 10, column = 1)
		ent_nwDiam.grid(row = 11, column = 1)
		btn_find_intersections.grid(row = 10, column = 2, rowspan = 2, sticky = tk.N + tk.S + tk.E + tk.W)

		#assigning variables to self so can be accessed by controller
		self.ent_initialTemp = ent_initialTemp
		self.ent_nwDiam = ent_nwDiam

		btn_show_graph = tk.Button(
			master = self,
			text = 'Show Graph',
			command = self.controller.set_plot
		)
		btn_show_graph.grid(row = 14, column = 0)

		btn_to_nanowire_mesh = tk.Button(
			master = self,
			text = 'Convert Graph to Nanowire Mesh',
			command = self.controller.to_nanowire_mesh
		)
		btn_to_nanowire_mesh.grid(row = 12, column = 2,
			rowspan = 2,
			sticky = tk.N + tk.S + tk.W + tk.E	
		)

		# scale bar info
		lbl_scale_bar_pixels = tk.Label(
			master = self,
			text = 'Scale Bar Pixels'
		)
		ent_scale_bar_pixels = tk.Entry(
			master = self
		)
		lbl_scale_bar_um = tk.Label(
			master = self,
			text = 'Scale Bar Length [um]'
		)
		ent_scale_bar_um = tk.Entry(
			master = self
		)

		lbl_scale_bar_pixels.grid(row = 12, column = 0)
		lbl_scale_bar_um.grid(row = 12, column = 1)
		ent_scale_bar_pixels.grid(row = 13, column = 0)
		ent_scale_bar_um.grid(row = 13, column = 1)
		#assigning variables to self so can be accessed by controller
		self.ent_scale_bar_pixels = ent_scale_bar_pixels
		self.ent_scale_bar_um = ent_scale_bar_um
		
		lbl_log = tk.Label(
			master = self,
			text = "Log"
		)
		lbl_log.grid(row = 14, column = 1)
		btn_print_parallel_edges = tk.Button(
			master = self,
			text = 'Print Parallel Edges',
			command = self.controller.print_parallel_edges
		)
		btn_print_parallel_edges.grid(row = 14, column = 2)
		self.txt_log = tk.Text(
			master = self
		)
		self.txt_log.grid(row = 15, column = 0, columnspan = 3, sticky = tk.N + tk.S + tk.E + tk.W)

