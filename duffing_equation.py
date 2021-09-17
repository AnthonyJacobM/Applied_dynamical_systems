import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.integrate import odeint
import PyDSTool as dst
from PyDSTool import *
from scipy.signal import argrelextrema
from scipy.fftpack import fft, fftn, ifft, ifftn

# change default matplotlib files
# dpi changes resolution of figures
mpl.rcParams['figure.dpi'] = 200
mpl.rcParams['savefig.dpi'] = 200
# fontsize is 18
mpl.rcParams['font.size'] = 18
# linewidth is 2
mpl.rcParams['lines.linewidth'] = 2.0

# references:
# https://en.wikipedia.org/wiki/Duffing_equation
# http://www.scholarpedia.org/article/Duffing_oscillator

# generate parameter dictionaries
# alpha ~ stiffness
# beta ~ nonlinearity of the restoring force
# delta ~ damping constant
# omega ~ angular frequency of forced periodic osicllation
# gamma ~ amplitude of the forced periodic oscillation

# traces
case_a = {'gamma1': 0.20, 'beta1': 1, 'delta1': 0.3, 'omega1': 1.2, 'alpha1': -1, 'x0': 1, 'v0': 0} # period 1 oscillation
case_b = {'gamma1': 0.28, 'beta1': 1, 'delta1': 0.3, 'omega1': 1.2, 'alpha1': -1, 'x0': 1, 'v0': 0} # period 2 oscillation
case_c = {'gamma1': 0.29, 'beta1': 1, 'delta1': 0.3, 'omega1': 1.2, 'alpha1': -1, 'x0': 1, 'v0': 0} # period 4 oscillation
case_d = {'gamma1': 0.37, 'beta1': 1, 'delta1': 0.3, 'omega1': 1.2, 'alpha1': -1, 'x0': 1, 'v0': 0} # period 5 oscillation
case_e = {'gamma1': 0.50, 'beta1': 1, 'delta1': 0.3, 'omega1': 1.2, 'alpha1': -1, 'x0': 1, 'v0': 0} # chaos


# cooler cases
case_cool_1 = {'gamma1': 8, 'beta1': 5, 'delta1': 0.02, 'omega1': 0.5, 'alpha1': 1, 'x0': 1, 'v0': 0} # chaos

def gen_sys(dict = case_a, tf = 100, N = 1000):
    """
    fucntion to generat the solution for the duffing equation!
    :param dict: dictionary containing parameter values ans initial conditions!
    :param tf: final time of simulation
    :param N: number of points for the numerical integration
    :return: z, t -- solution of x, y, and time integrated
    """
    # unpack the parameters
    gamma1 = dict['gamma1']
    beta1 = dict['beta1']
    alpha1 = dict['alpha1']
    omega1 = dict['omega1']
    delta1 = dict['delta1']

    # initial conditions for the simulation
    x0 = np.array([dict['x0'], dict['v0']])

    # generate a time sequence
    t = np.linspace(0, tf, N)

    # generate the ode
    def sys(X, t = 0):
        """
        ode system of the duffing equation!
        :param X: array data containing two variables -- position and velocity
        :param t: time
        :return: Y -- integrated solution
        """
        # generate an array for the integration data
        Y = np.zeros(len(X))
        # notation:
        # X[0] ~ x, X[1] ~ y
        # Y[0] ~ x', Y[1] ~ y'
        Y[0] = X[1]
        Y[1] = gamma1 * np.cos(omega1 * t) - (delta1 * X[1] + X[0] * (alpha1 + beta1 * X[0]**2))

        return Y

    # integrate the solution....
    z, infodict = odeint(sys, x0, t, full_output = True) # integrated solution

    return z, t

def plot_solution(dict = case_a, tf = 100, N = 1000, xy_bool = True, t_bool = True, dmap_bool = True, fourier_bool = True):
    """
    function to plot the time series solution of the duffing equation
    :param dict: dicitionary of the cases
    :param tf: final time of simulation
    :param N: Number of points used for numerical integration
    :param xy_bool: boolean for the plot of position versus velocity
    :param t_bool: boolean for the time plot
    :return: plotted solution
    """
    # genrate a figure
    # integrate the solution
    z, t = gen_sys(dict = dict, tf = tf, N = N)

    # unpack key variables
    x, v = z.T

    # plot
    if xy_bool == True:
        fig = plt.figure()
        # generate quvier data!
        dist = 0.15 # distance of the quiver plot
        Nquiv = 15 # number of quiver arrows
        xlow = (1 - dist) * np.min(x) # minimum of the data
        xhigh = (1 + dist) * np.max(x) # maximum of the data
        vlow = (1 - dist) * np.min(v) #...
        vhigh = (1 + dist) * np.max(v) #...
        x1= np.linspace(xlow, xhigh, Nquiv)
        v1 = np.linspace(vlow, vhigh, Nquiv)

        # generate a meshgrid
        x1mesh, v1mesh = np.meshgrid(x1, v1)

        # generate arrow over the vector field
        # unpack the parameters
        gamma1 = dict['gamma1']
        beta1 = dict['beta1']
        alpha1 = dict['alpha1']
        omega1 = dict['omega1']
        delta1 = dict['delta1']
        def gen_arrows(X, t=0):
            """
            ode system of the duffing equation!
            :param X: array data containing two variables -- position and velocity
            :param t: time
            :return: Y -- integrated solution
            """
            # generate an array for the integration data
            # notation:
            # X[0] ~ x, X[1] ~ y
            # Y[0] ~ x', Y[1] ~ y'
            x1 = X[1]
            x2 = gamma1 * np.cos(omega1 * t) - (delta1 * X[1] + X[0] * (alpha1 + beta1 * X[0] ** 2))

            return np.array([x1, x2])
        dx1, dv1 = gen_arrows([x1mesh, v1mesh])

        # hypot associated to the data arrows!
        M = (np.hypot(dx1, dv1))
        M[M == 0] = 1 # no division by zeros!
        dx1 /= M # normalize the arrows
        dv1 /= M # ...

        # generate the quiver plot
        plt.quiver(x1, v1, dx1, dv1, M, pivot = 'mid', cmap = 'RdBu')
        plt.plot(x,v, 'k-', lw = 2)
        plt.xlabel('Position', fontsize = 18)
        plt.ylabel('Velocity', fontsize = 18)
        plt.show()

    if t_bool == True:
        fig = plt.figure()
        plt.plot(t, x, 'b-', lw = 2)
        plt.xlabel('t', fontsize = 18)
        plt.ylabel('Position', fontsize = 18)
        plt.show()

    if dmap_bool == True:
        idx_peaks = argrelextrema(x, np.greater)[0] # indices of the peaks
        idx_valleys = argrelextrema(x, np.less)[0] # indices of the valleys
        peaks = x[idx_peaks] # peaks of the position
        valleys = x[idx_valleys] # troughs of the position
        n_peaks = len(peaks)
        n_valleys = len(valleys)
        peaks_x = np.linspace(0, n_peaks, n_peaks)
        valleys_x = np.linspace(0, n_valleys, n_valleys)

        # plot it!
        plt.plot(peaks_x, peaks, 'k', lw = 2)
        plt.xlabel('n', fontsize = 18)
        plt.ylabel('Peaks', fontsize = 18)
        plt.show()

        plt.plot(valleys_x, valleys, 'k', lw = 2)
        plt.xlabel('n', fontsize = 18)
        plt.xlabel('Valleys', fontsize = 18)
        plt.show()

    if fourier_bool == True:
        # generate fourier coefficients
        z = fft(x)
        # power
        mag = np.abs(z)

        plt.plot(mag, 'k', lw = 2)
        plt.xlabel('f', fontsize = 2)
        plt.ylabel('Magnitude', fontsize = 18)
        plt.show()



def run_simulation(tf = 100, N = 1000):
    """
    function to run simulation
    :return:
    """
    flag = 1
    while flag != 0:
        option = input('choose: case_a, case_b, case_c, case_d, case_e, case_cool_1, otherwise 0 to quit')
        if option == 'case_a':
            dict = case_a
        elif option == 'case_b':
            dict = case_b
        elif option == 'case_c':
            dict = case_c
        elif option == 'case_d':
            dict = case_d
        elif option == 'case_e':
            dict= case_e
        elif option == 'case_cool_1':
            dict = case_cool_1
        else:
            flag = 0
        plot_solution(dict = dict, tf = tf, N = N)



run_simulation()