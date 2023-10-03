from matplotlib import pyplot as plt
def _idx_to_xy(idx):
	x,y = idx[1], idx[0]
	return x,y

def _xy_to_idx(xy):
	row, col = xy[1], xy[0]
	return row,col

def plot_imgs(imgs, colorbar = False):
	fig, axs = plt.subplots(
		nrows = (rows := round(len(imgs)/2 + 0.1)),
		ncols = 2,
		figsize = (14, 7 * rows)
	)
	for n, (title, img) in enumerate(imgs.items()):
		ax = axs.flatten()[n]
		ax.set_title('{}'.format(title))
		plot = ax.imshow(
			img,
			cmap = 'hot'
		)
		if colorbar:
			cbar = fig.colorbar(
				plot,
				ax = ax
			)
			cbar.minorticks_on()

	plt.show()
