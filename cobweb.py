# import relevant libraries
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

import scipy.signal
from matplotlib import rc
from mkfolder import create_folder
# Use LaTeX throughout the figure for consistency
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 16})
rc('text', usetex=True)
# Figure dpi
dpi = 72

path_base = r'C:\Users\antho\Documents\Classes\Fall 2021\Applied Dynamical Systems\Python'

#create_folder('logistic')
#create_folder('time_evolution')
#create_folder('cweb')

def plot_cobweb(f, r, x0, nmax = 40, x_low = 0, x_high = 1, n_x = 500, save_bool = False, show_bool = True, model = 'logistic'):
    """
    function to plot a cobweb
    :param f: function
    :param r: parameter used for the function
    :param x0: initial condition
    :param nmax: max number of cobweb iterations
    :param x_low: minimum on the x array
    :param x_high: maximum on the x array
    :param n_x: length of the x array
    :param show_bool: boolean to show the figure
    :param save_bool: boolean to save the current figure
    :param model: model used for the time series
    :return: plotted function
    """
    # generate linear spaced points
    x = np.linspace(x_low, x_high, n_x)
    # generate a figure
    plt.close()
    plt.clf()
    fig = plt.figure()

    # plot the function
    plt.plot(x, f(x, r), color = 'black', lw = 2)
    # plot the 45 degree line
    plt.plot(x, x, color = 'black', lw = 2)

    # perform an iteration
    px, py = np.empty((2, nmax + 1, 2))
    # set the initial values of the iteration
    px[0], py[0] = x0, 0
    # get all the iterations over the discrete map
    for n in range(1, nmax, 2):
        px[n] = px[n-1]
        py[n] = f(px[n-1], r)
        px[n+1] = py[n]
        py[n+1] = py[n]

    # plot the trace over the iteration
    plt.plot(px, py, color = 'b', lw = 1, alpha = 0.5)
    plt.grid(which='minor', alpha=0.5)
    plt.grid(which='major', alpha=0.5)
    plt.xlabel(r'$x_n$', fontsize = 18)
    plt.ylabel(f.latex_label, fontsize = 18)
    plt.title(r'$x_0 = {}, r = {}$'.format(np.round(x0, 3), np.round(r, 3)))
    plt.tight_layout()
    if show_bool == True:
        plt.show()


    if save_bool == True:
        file_name = 'cobweb_{}_{}'.format(np.round(x0*10000, 3), np.round(r*10000, 3))
        if model == 'logistic':
            path = r'C:\Users\antho\Documents\Classes\Fall 2021\Applied Dynamical Systems\Python\figures\logistic\cweb'
        plt.savefig(path + f"\{file_name}.png", dpi = 300)

    # plot the trace over the iteration
    plt.close()
    plt.clf()
    fig = plt.figure()
    plt.plot(py, color = 'b', lw = 1, alpha = 0.5)
    plt.grid(which='minor', alpha=0.5)
    plt.grid(which='major', alpha=0.5)
    plt.xlabel(r'$n$', fontsize = 18)
    plt.ylabel(f.latex_label, fontsize = 18)
    plt.title(r'$x_0 = {}, r = {}$'.format(np.round(x0, 3), np.round(r, 3)))
    plt.tight_layout()
    if show_bool == True:
        plt.show()

    if save_bool == True:
        file_name = 'time_{}_{}'.format(np.round(x0*10000, 3), np.round(r*10000, 3)).strip('.')
        if model == 'logistic':
            path = r'C:\Users\antho\Documents\Classes\Fall 2021\Applied Dynamical Systems\Python\figures\logistic\time_evolution'
        plt.savefig(path + f"\{file_name}.png", dpi = 300)

    return py

class AnnotatedFunction:
    """A small class representing a mathematical function.
    This class is callable so it acts like a Python function, but it also
    defines a string giving its latex representation.
    """

    def __init__(self, func, latex_label):
        self.func = func
        self.latex_label = latex_label

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

# The logistic map, f(x) = rx(1-x).
func = AnnotatedFunction(lambda x,r: r*x*(1-x), r'$rx(1-x)$')
#plot_cobweb(func, 2.8, 0.2, save_bool = True)
#plot_cobweb(func, 3.8, 0.2, 200, save_bool = True)

#create_folder(name = folder_name)
def generate_cobweb_figures(f, r, x0, nmax, x_low, x_high, n_x, r_low, r_high, n_r, model = 'logistic'):
    """
    function to generate a movie
    :param f: function used for the iterative map
    :param r: parameter
    :param x0: initial condition
    :param nmax: maximum number of iterations
    :param x_low: lower  bound on x
    :param x_high: upper bound on x
    :param n_x: number of x elelements in x
    :param r_low: lower bound on r
    :param r_high: upper bound on r
    :param n_r: number of iterations over r
    :param folder_name: folder used to generate a sequence of plots
    :return: collection of plotted figures
    """
    # create an array over the parameter
    r_array = np.linspace(r_low, r_high, n_r)
    # iterate over the sequence
    for i, i0 in enumerate(r_array):
        i0 = np.round(i0, 4)
        print('r = ' + str(i0))
        rval = i0
        plot_cobweb(f, rval, x0, nmax=40, x_low=0, x_high=1, n_x=500, save_bool=True, show_bool=False, model = model)



#generate_cobweb_figures(f = func, r = 0, x0 = 0.2, nmax = 40, x_low = 0, x_high = 1, n_x = 500, r_low = 0, r_high = 4, n_r = 100, model = 'logistic')


import moviepy
import moviepy.video.io.ImageSequenceClip
def make_movie(fps = 15, path = path_base, output_name = 'my_movie', type = 'cobweb', model = 'logistic'):
    """
    function to make a movie of pertrubed time series about the bifurcation of the parameter: par
    :param fps: frames per second: 30, 60, etc.
    :param image_folder: folder to save the images as save .png files...
    :param output_name: file name of the saved movie.
    :param type: 'time_series', 'reff', 'ib_sg', 'v_ib', 'ib_sb', 'sb_sg', 'bifurcation'
    :return: movie (saved)
    """
    if type == 'cobweb':
        type = 'cweb'
    elif type == 'time':
        type = 'time_evolution'

    image_folder = r'C:\Users\antho\Documents\Classes\Fall 2021\Applied Dynamical Systems\Python\figures\logistic' + f"\{type}"
    if output_name == 'my_movie':
        output_name = f"{type}_{model}_movie"

    image_files = [image_folder + f"\{img}" for img in os.listdir(image_folder) if img.endswith("png")]
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(f"{output_name}.mp4")





#make_movie(type = 'time')



# genererate a curve for a linear combination of curves....
x = np.linspace(0, 10, 101)
y1 = 2 * np.sin(np.pi * x * 3 - 1)
y2 = -1.5 * np.cos(np.pi * x * 2 - 3)
y3 = 3 * np.sin(np.pi * x * 5 + 2)
y = y1 + y2 + y3
sine_curve = AnnotatedFunction(lambda x: (2 * np.sin(np.pi * x * 3 - 1) -1.5 * np.cos(np.pi * x * 2 - 3) + 3 * np.sin(np.pi * x * 5 + 2)) / 6, r'Sine Curve')

plt.plot(x, y, 'k-', lw = 2)
plt.show()




# generate the cobweb figure associated to the time series of the data.....
def plot_cobweb_sine(sine_curve, x0, nmax = 40, x_low = 0, x_high = 10, n_x = 1000, save_bool = False, show_bool = True):
    """
    function to plot a cobweb
    :param f: function
    :param r: parameter used for the function
    :param x0: initial condition
    :param nmax: max number of cobweb iterations
    :param x_low: minimum on the x array
    :param x_high: maximum on the x array
    :param n_x: length of the x array
    :param show_bool: boolean to show the figure
    :param save_bool: boolean to save the current figure
    :return: plotted function
    """
    # generate linear spaced points
    x = np.linspace(x_low, x_high, n_x)
    # generate a figure
    plt.close()
    plt.clf()
    fig = plt.figure()

    # plot the function
    plt.plot(x, sine_curve(x), color = 'black', lw = 2)
    # plot the 45 degree line
    plt.plot(x, x, color = 'black', lw = 2)

    # perform an iteration
    px, py = np.empty((2, nmax + 1, 2))
    # set the initial values of the iteration
    px[0], py[0] = x0, 0
    # get all the iterations over the discrete map
    for n in range(1, nmax, 2):
        px[n] = px[n-1]
        py[n] = sine_curve(px[n-1])
        px[n+1] = py[n]
        py[n+1] = py[n]

    # plot the trace over the iteration
    plt.plot(px, py, color = 'b', lw = 1, alpha = 0.5)
    plt.grid(which='minor', alpha=0.5)
    plt.grid(which='major', alpha=0.5)
    plt.xlabel(r'$x_n$', fontsize = 18)
    plt.ylabel(sine_curve.latex_label, fontsize = 18)
    plt.title(r'$x_0 = {}$'.format(np.round(x0, 3)))
    plt.tight_layout()
    if show_bool == True:
        plt.show()



    # plot the trace over the iteration
    plt.close()
    plt.clf()
    fig = plt.figure()
    plt.plot(py, color = 'b', lw = 1, alpha = 0.5)
    plt.grid(which='minor', alpha=0.5)
    plt.grid(which='major', alpha=0.5)
    plt.xlabel(r'$n$', fontsize = 18)
    plt.ylabel(sine_curve.latex_label, fontsize = 18)
    plt.tight_layout()
    if show_bool == True:
        plt.show()


    return py


plot_cobweb_sine(sine_curve, x0 = 0, x_high = 1, x_low = -1)