#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Showing the importance of the rate determining step. The global speed is shown as well as the speed of the rate determinig step.o

The mecanism proposed is 

A ->k1 B ->k2 C ->k3 D, 
    all reactions being first order
    Initial condition : only A is present
    k1 = k3 = 100
    k2 is varied according to the slider.

Informations
------------
Author : Martin Vérot  from the ENS de Lyon, France, with the help of some scripts to create the buttons and sliders taken from https://github.com/araoux/python_agregation (written by P Cladé, A Raoux and F Levrier)
Licence : Creative Commons CC-BY-NC-SA 4.0

WARNING this program requires the widgets.py file to work
"""

# Importation of libraries
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.integrate
import widgets
#plt.rc('text', usetex=True)
#plt.rc('text.latex', preamble=r'\usepackage{stackrel}')

# Definition of functions 
def kineticEquations(Concentrations, t, k1,k2,k3):
    """
        Concentrations is an array containing the concentrations of A, B, C, D in that order
        t is the time
        k1,k2,k3 are the kinetic constants
    """
    #print(Concentrations)
    KineticMatrix = np.zeros((Concentrations.shape[0],Concentrations.shape[0]))
    ReactionMatrix = np.array([[-1,1,0,0],[0,-1,1,0],[0,0,-1,1],[0,0,0,0]]).transpose()
    KineticConstants = np.array([k1,k2,k3,0.])
    KineticMatrix = ReactionMatrix*KineticConstants
    finalMatrix =  np.matmul(KineticMatrix,Concentrations) 
    #if t<0.03:
    #    print(KineticConstants)
    #    print(ReactionMatrix)
    #    print(KineticMatrix)
    #    print('final')
    #    print(finalMatrix)
    #print(finalMatrix)
    return(finalMatrix)


def plot_data(ratio):
    k2 = k1*10**(ratio)
    params=(k1,k2,k3)
    t = np.linspace(0,10/min(params),10000) 
    vals = scipy.integrate.odeint(kineticEquations,C0,t ,args = params)
    Concentrations = vals.transpose()
    v2 = k2*Concentrations[1]
    v3 = k3*Concentrations[2]
    ax1.set_xlim(min(t),max(t))
    ax2.set_xlim(min(t),max(t))
    ax3.set_xlim(min(t),max(t))
    ax4.set_xlim(0,max(v2)*1.05)
    ax1.set_ylim(0,max(max(Concentrations[0]),max(Concentrations[3])))
    ax2.set_ylim(0,max(max(Concentrations[1]),max(Concentrations[2])))
    ax3.set_ylim(0,max(max(k2*Concentrations[1]),max(k3*Concentrations[2])))
    ax4.set_ylim(0,1.05*max(v3))
    
    lines['A1'].set_data(t,Concentrations[0])
    lines['B1'].set_data(t,Concentrations[1])
    lines['C1'].set_data(t,Concentrations[2])
    lines['D1'].set_data(t,Concentrations[3])
    lines['A2'].set_data(t,Concentrations[0])
    lines['B2'].set_data(t,Concentrations[1])
    lines['C2'].set_data(t,Concentrations[2])
    lines['D2'].set_data(t,Concentrations[3])
    lines['vd'].set_data(t,v3)
    lines['vb'].set_data(t,v2)
    lines['versus'].set_data(v2,v3)
    lines['straight'].set_data(v2,v2)
    fig.canvas.draw_idle()

# Main program
if __name__ == "__main__":
    #kinetic constants
    k1=100 #s^-1
    k2=1   #s^-1
    k3=100 #s^-1
    #Initial conditions : only A is present
    C0 = [1.,0.,0.,0.]



    parameters = {
        'ratio' : widgets.FloatSlider(value=-2, description='$\log(k_2/k_1)$', min=-3, max=3),
    }
    fig,axes=plt.subplots(2,2)
    ax1 = plt.subplot(2,2,1)
    ax2 = plt.subplot(2,2,2)
    ax3 = plt.subplot(2,2,3)
    ax4 = plt.subplot(2,2,4)
    ax1.set_title('Concentrations versus time')
    ax2.set_title('Concentrations versus time')
    #ax2.set_title('B')
    ax3.set_title('$v_d(\mathrm{B})$ and $v_a(\mathrm{D})$ versus time')
    ax4.set_title('$v_a(\mathrm{D})=f(v_d(\mathrm{B}))$')
    #ax4.set_title('D')

    lines = {}
    lines['A1'], = ax1.plot([],[],label='[A]',color='C0')
    lines['B1'], = ax1.plot([],[],label='[B]',color='C1',ls='--')
    lines['C1'], = ax1.plot([],[],label='[C]',color='C2',ls='--')
    lines['D1'], = ax1.plot([],[],label='[D]',color='C3')
    lines['A2'], = ax2.plot([],[],label='[A]',color='C0',ls='--')
    lines['B2'], = ax2.plot([],[],label='[B]',color='C1')
    lines['C2'], = ax2.plot([],[],label='[C]',color='C2')
    lines['D2'], = ax2.plot([],[],label='[D]',color='C3',ls='--')
    lines['vd'], = ax3.plot([],[],label='v_a(D)',color='C3')
    lines['vb'], = ax3.plot([],[],label='v_d(B)',color='C2')
    lines['versus'], = ax4.plot([],[],color='C4')
    lines['straight'], = ax4.plot([],[],color='C5')
    ax1.legend(loc='right')
    ax2.legend(loc='right')

    fig.suptitle(r'Rate determining step for a reaction $A \rightarrow^{k_1} B \rightarrow^{k_2} C \rightarrow^{k_3} D$ with $k_1=k_3$', fontsize=16)
    param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.35, 0.91, 0.4, 0.02])



    plt.show()
