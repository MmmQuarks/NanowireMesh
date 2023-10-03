import tkinter as tk
from Model import Model #importing with relative path
from ImageView import ImageView
from ActionView import ActionView
from Controller import Controller

class Window(tk.Tk):
	def __init__(
		self, 
		labeled_img = None,
		profiles_and_img = None,
		graph = None,
		**kwargs
	):
		super().__init__(**kwargs)

		self.title('Nanowire Image Reconstruction')

		# create model
		model = Model()
		self.model = model
  
		# creating container frames for our custom objects
		l_frame, r_frame = [tk.Frame(master = self) for n in range(2)]
		l_frame.grid(row = 0, column = 0, sticky = 'nsew', padx = 20)
		r_frame.grid(row = 0, column = 1, sticky = 'nsew')
		self.rowconfigure(0, weight = 1, minsize = 500 )
		self.columnconfigure(0, weight = 9,  minsize = 400)
		self.columnconfigure(1, weight = 1,  minsize = 50)

		# creating controller
		controller = Controller(model)
		self.controller = controller

		# creating an image view 
		image_view = ImageView(master = l_frame, controller = controller)
		# letting the controller know it now has an ImageView
		controller.set_image_view(image_view)
		# adding the image view to the left frame
		image_view.pack(side = tk.TOP, fill = tk.BOTH, expand = True)

		self.image_view = image_view

		# creating an action view
		action_view = ActionView(master = r_frame, controller = controller)
		# letting the controller know it now has an ImageView
		controller.set_action_view(action_view)
		# adding the action_view to the right frame
		action_view.pack(side = tk.TOP, fill = tk.BOTH, expand = True)

		if profiles_and_img is not None:
			controller.load_profiles_and_img(file_path = profiles_and_img)
		if labeled_img is not None:
			controller.load_labeled_img(file_path = labeled_img)
		if graph is not None:
			controller.load_graph(file_path = graph)