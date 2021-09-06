import PyDSTool as dst
from PyDSTool import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# change default matplotlib files
# dpi changes resolution of figures
mpl.rcParams['figure.dpi'] = 200
mpl.rcParams['savefig.dpi'] = 200
# fontsize is 18
mpl.rcParams['font.size'] = 18
# linewidth is 2
mpl.rcParams['lines.linewidth'] = 2.0

# reference:
# https://en.wikipedia.org/wiki/Pitchfork_bifurcation

def generate_system(tf = 100, case = 'subcritical', integrator = 'Vode'):
    """
    function to generate the pitchfork bifurcation
    :param tf: final time of integration
    :param case: subcritical or supercritical
    :param integrator: Radau or Vode
    :return: ode -- generated system
    """
    # initialize the integrator
    if integrator == 'Radau':
        gentype = dst.Generator.Radau_ODEsystem # generator
    else:
        gentype = dst.Generator.Vode_ODEsystem # generator

    DSargs = dst.args(name='pitchfork')  # generate the system
    if case == 'subcritical':
        sgnvalue = 1
    else:
        sgnvalue = -1
    DSargs.pars = {'r': 0, 'sgn': sgnvalue}  # parameters
    DSargs.varspecs = {'x': 'r * x + sgn * x*x*x'}  # rhs of the differential equations
    DSargs.ics = {'x': 1e-3}  # initial conditions
    DSargs.tdomain = [0, tf]  # time domain

    # generate the system
    ode = gentype(DSargs)

    return ode

def generate_bifurcation(ode, show_bool = True, title_bool = True):
    """
    function to generate the bifurcation for the pitchfork bifurcation
    :param ode: ode -- generated from generate_system()....
    :param show_bool: boolean to plot the bifurcation
    :param title_bool: boolean to show the title of the figure
    :return: PC -- python continuation curve for the system
    """
    # generate a bifurcation
    PC = ContClass(ode)  # generate a continuation class...
    PCargs = dst.args(name='EQ', type='EP-C')
    PCargs.MaxStepSize = 0.01
    PCargs.MinStepSize = 1e-6
    PCargs.freepars = ['r']
    PCargs.LocBifPoints = 'all'
    PCargs.SaveEigen = True
    PC.newCurve(PCargs)

    # begin the continuation:
    PC['EQ'].forward()
    PC['EQ'].backward()

    PCargs = dst.args(name='EQnew', type='EP-C')
    PCargs.initpoint = {'x': 0, 'r': 1}
    PCargs.MaxStepSize = 0.01
    PCargs.MinStepSize = 1e-6
    PCargs.freepars = ['r']
    PCargs.LocBifPoints = 'all'
    PCargs.SaveEigen = True
    PC.newCurve(PCargs)
    PC['EQnew'].forward()
    PC['EQnew'].backward()

    if show_bool == True:
        PC.display(['r', 'x'], stability=True)
        if title_bool == True:
            par_dict = ode.pars
            if par_dict['sgn'] > 0:
                title = 'Subcritical'
            else:
                title = 'Supercritical'

            plt.title(title)
        plt.show()

    return PC

def plot_time_series(ode, tf = 100, user_input_bool = True, title_bool = True):
    """
    function to plot the time series of the pitchfork bifurcation
    :param ode: ode -- returned from generate_system()
    :param tf: final time of integration
    :param user_input_bool: boolean to request user inputed intial condition and parameter value
    :param title_bool: boolean to plot with the title
    :return: plotted solution in mpl figure
    """
    flag = 1
    pts = ode.compute('traj').sample(dt = 0.1)
    pardict = ode.pars
    r0 = pardict['r']
    if r0 > 0:
        col = 'red'
    else:
        col = 'blue'
    plt.plot(pts['t'], pts['x'], color = col)
    plt.xlabel('t')
    plt.ylabel('x')
    while flag != 0:
        if user_input_bool == True:
            continue_bool = int(input('Type "1" to continue and "0" to stop'))
            flag = int(continue_bool)
            print(continue_bool)
            if continue_bool == True or continue_bool == 1:
                    r0 = float(input('choose value of r'))
                    x0 = float(input('choose value of x'))
                    ode.set(pars = {'r': r0})
                    ode.set(ics = {'x': x0})
                    pts = ode.compute('test' + str(r0) + str(x0)).sample(dt = 0.1)
                    if r0 > 0:
                        col = 'red'
                    else:
                        col = 'blue'
                    plt.plot(pts['t'], pts['x'], 'k-', color = col)
            else:
                break
        else:
            break

    if title_bool == True:
        par_dict = ode.pars
        if par_dict['sgn'] > 0:
            title = 'Subcritical'
        else:
            title = 'Supercritical'

        plt.title(title)

    plt.show()


def run_simulation(case = 'subcritical', tf = 35):
    """
    function to run the simulation
    :param case: "supercritical" or "subcritical"
    :param tf: final time of the simulation
    :return: plotted bifurcations and time series
    """
    if case == 'subcritical':
        ode_sub = generate_system(tf=tf, case='subcritical')
        PC_sub = generate_bifurcation(ode_sub)
        plot_time_series(ode = ode_sub)

    else:
        ode_super = generate_system(tf=tf, case='supercritical')
        PC_sup = generate_bifurcation(ode_super)
        plot_time_series(ode=ode_super)

run_simulation(case = 'subcritical')
run_simulation(case = 'supercritical')


