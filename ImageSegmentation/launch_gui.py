#import pickle
#import tkinter as tk
#import skimage
#import numpy as np
#from PIL import Image, ImageTk
#from matplotlib.lines import Line2D
#import NanowireImageSegmentation as nim
#from matplotlib import pyplot as plt
#from LabeledImage import LabeledImage

import sys
import pdb
from Window import Window
import argparse






def main():
	if len(sys.argv) > 1:
		parser = argparse.ArgumentParser()
		group = parser.add_mutually_exclusive_group()
		group.add_argument(
			'--labeled-img',
			help = 'path to labeled image to load',
			default = None
		)
		group.add_argument(
			'--profiles-and-img',
			help = 'path to pickle of profiles and image to load. profiled image will be clustered upon opening',
			default = None
		)
		parser.add_argument(
			'--graph',
			help = 'path to graph to load',
			default = None
		)
		args = parser.parse_args()
		window = Window(
			labeled_img = args.labeled_img,
			profiles_and_img = args.profiles_and_img,
			graph = args.graph
		)
	else:
		window = Window()

	window.mainloop()
	pdb.set_trace()

if __name__ == '__main__':
	main()