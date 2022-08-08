#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : 2020-2023
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
Modified by   : Jeremy Lapierre
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""

import tempfile
import os
import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import argparse

class CTE:
    """Some importnat physical constants
    """
    Kb = 8.314462618E-3 #kJ/(mol⋅K) (kNA)

def KbT(absolute_temperature):
    """Return the value of Kb*T in kJ/mol

    Args:
        absolute_temperature (float): The absolute temperature in kelvin
    """
    
    return absolute_temperature*CTE.Kb

def beta(absolute_temperature):
    return 1 / KbT(absolute_temperature)

def run(command, shell = True, executable = '/bin/bash', Popen = False):
    #Here I could make some modification in order that detect the operator system
    #NAd make the command compatible with the opertor system
    #the fucntion eval could be an option if some modifcation to the variable command 
    #need to be done.... SOme fligth ideas...

    if Popen:
        #In this case you could acces the pid as: run.pid
        process = subprocess.Popen(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    else:
        process = subprocess.run(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
        returncode = process.returncode
        if returncode != 0:
            print(f'Command {command} returned non-zero exit status {returncode}')
            raise RuntimeError(process.stderr)
    return process


def PandeWeights(temperatures, energies, alpha = None, weight_GROMACS_format = True):
    # Here explain the definition of the gromacs weight
    # https://manual.gromacs.org/documentation/2019/reference-manual/algorithms/expanded-ensemble.html
    # The use of the alpha is not completly right from a theretical point of view for this method
    assert len(temperatures) == len(energies), "temperatures and energies must have the same number of elements"
    # Deffining the coefficient matrix
    M = 1*np.eye(len(temperatures), k = 1, dtype=int)
    np.fill_diagonal(M, -1)
    M[-1] = np.ones(len(temperatures), dtype = int)
    # Solution Vector
    b = np.zeros(len(temperatures))
    # Initialize b[0] dependign if alpha was selected
    if alpha:
        gamma = alpha/(1-alpha)
        b[0] = 0.5*(beta(temperatures[1]) - beta(temperatures[0]))*(energies[0] + energies[1]) - np.log(gamma*(len(temperatures) - 1))
    else:
        b[0] = 0.5*(beta(temperatures[1]) - beta(temperatures[0]))*(energies[0] + energies[1])

    for i in range(1, len(energies)):
        try:
            b[i] = 0.5*(beta(temperatures[i + 1]) - beta(temperatures[i]))*(energies[i] + energies[i + 1])
        # I reached the last element of the list so I need to add cero to the array but I already initialized on zero.
        except:
            break
    # Solve the system of equations
    # Calculate the weights:
    weights = np.linalg.solve(M, b)

    # Convert if needed to GROMACS weights
    if weight_GROMACS_format:
        minimum = weights.min()
        if minimum > 0:
            weights -= minimum
        else:
            weights += np.abs(minimum)
    return weights

def get_energy(edr, annealing_times, energy_type = 'Potential'): # Could be Total-Energy
    """The two input is a list that could be obtained from mdp.annealing3
    The structure is the following
    annealing_times = [t1, t2, t3, t4, t5, t6]
    These mean:
    From mdp.annealing, we also obtain the list
    annealing_temperatures = [T1, T1, T2, T2, T3, T3]
    Tis two list are connected and mean:
    from time t1 to t2 the temperature was constant, therefore from here we should get the average energy for T1.
    Then from t2 to t3 the temperature increase from T1 to T2
    Then from t3 to t4 the temperature is constant at T2,therefore from here we should get the average energy for T2.
    And so on and so for.

    Args:
        annealing_temperatures (_type_): _description_
        annealing_times (_type_): _description_
    """
    fig, ax = plt.subplots(figsize = (16,9))
    data = pd.DataFrame()
    xvg_tmp_file = tempfile.NamedTemporaryFile(suffix='.xvg')
    energy = []
    iterator = range(0, len(annealing_times)-1, 2)

    for state, index in tqdm.tqdm(enumerate(iterator), total=len(iterator)):#enumerate(iterator):# # the calculation is per pair of times, beetween the first to time the temperature was keep constant, then the system was heated and repeated again.
        result = run(f"export GMX_MAXBACKUP=-1; echo {energy_type} | gmx energy -f {edr} -b {annealing_times[index]} -e {annealing_times[index + 1]} -o {xvg_tmp_file.name} | grep \'{energy_type.replace('-',' ')}\'")
        energy.append(float(result.stdout.split()[-5]))
        """
        Energy                      Average   Err.Est.       RMSD  Tot-Drift
        -------------------------------------------------------------------------------
        Potential                -1.30028e+06         --     1682.1   -2422.24  (kJ/mol)
        Total Energy                -952595         --    2606.81    -3688.3  (kJ/mol)
        """

    return energy

def get_weights_from_log(log, plot = False):
    """This function read a log file obtained from simulated tempering and returns a numpy array
    of shape (number of time points, number of states, 2)
    The first axis of the array refears to each output on the log file.
    The second axis refears to the temperature states
    And the last one to the counts and value of the corresponded weights
    E.g., if we would like to know the value of the weight on the second output of the log file for the temperature state 5.
    >>> from mdynamic.simutemp import weights
    >>> time, info_weights = weights.get_weights_from_log('tempering.log')
    >>> info_weights[1,4,1]
    387.69299

    Get the final weights:
    >>> info_weights[-1,:,1]
    array([   0.     ,   95.78638,  194.94769,  291.73962,  387.69299,
        484.79767,  563.36194,  609.72424,  656.0224 ,  702.54724,
        747.32648,  793.23291,  838.1333 ,  880.62061,  922.65308,
        966.77673, 1008.16711, 1051.17847, 1091.65161, 1133.02808,
       1175.10168, 1212.64014, 1249.82324, 1289.07568, 1329.84229,
       1366.30273])

    Args:
        log (str): path to the log file

    Returns:
        numpy.array: The weights information.
    """
    with open(log, 'r') as f:
        log_file = f.readlines()

    i = 0
    time = []
    weights_info = []
    weights_0 = []
    while i < len(log_file):
        if 'init-lambda-weights[' in log_file[i]:
            weights_0.append(float(log_file[i].split('=')[-1]))

        if 'MC-lambda information' in log_file[i]:
            # Finding the time
            for j in range(i,0,-1):
                if log_file[j].startswith('           Step           Time'):
                    j += 1
                    time.append(float(log_file[j].split()[-1]))
                    break
            # Finding the weight
            weights_info_tmp = []
            i += 3
            while log_file[i] != '\n':
                split = log_file[i].split()
                count = int(split[2])
                weight = float(split[3])
                weights_info_tmp.append((count, weight))
                i += 1
            weights_info.append(weights_info_tmp)
        i += 1
    # Add weights at t = 0, because the counts are all 0 and I delate the entrances with total count 0 in next lines,
    # What i could do is put 1 in the initial temperature
    time.insert(0,0)
    weights_info.insert(0,list(zip([1] + (len(weights_0) - 1)*[0], weights_0)))

    #Converting to array
    time = np.array(time)
    weights_info = np.array(weights_info)
    # Some times (I don't know why) GROMACS reset all the weights and all the counts are 0. We need to eliminate those points
    sum_of_weights = weights_info[:,:,0].sum(axis = 1)
    time = time[sum_of_weights != 0]
    weights_info = weights_info[sum_of_weights != 0]
    sum_of_weights = sum_of_weights[sum_of_weights != 0]


    if plot:
        dir = os.path.dirname(log)
        fig, axes = plt.subplots(2, figsize = (16,9), sharex=True)
        NUM_COLORS = weights_info.shape[1]
        cm = plt.get_cmap('viridis')#gist_rainbow viridis
        for axe in axes:
            axe.set_prop_cycle('color', [cm(1.*j/NUM_COLORS) for j in range(NUM_COLORS)])

        probability = weights_info[:,:,0] / sum_of_weights[:,np.newaxis]
        for j in range(weights_info.shape[1]):
            #axes[0].plot(time, weights_info[:,j,0], label = str(j))
            axes[0].plot(time, probability[:,j], label = str(j))
            axes[1].plot(time, weights_info[:,j,1])
        
        fig.legend(loc = 'lower center', ncol = int(weights_info.shape[1] / 2))
        axes[0].set(
            xlim = (time.min(), time.max()),
            ylim = (0,1),
            ylabel = 'Probability',
            )
        axes[1].set(
            xlabel = 'Time [ps]',
            ylabel = 'Weight values'
        )
        #plt.show()
        fig.savefig(os.path.join(dir,'weights_progression.svg'), bbox_inches="tight")

        ## Plotting the violin plot of the weights
        #df = pd.DataFrame()
        #for j in range(weights_info.shape[1]):
        #    #df[temperatures[j]] = weights_info[:,j,1]
        #    df[j] = weights_info[:,j,1]
        ## Set up the matplotlib figure
        #sns.set_theme(style="whitegrid")
        #fig, ax = plt.subplots(figsize=(25, 25))

        ## Draw a violinplot with a narrower bandwidth than the default
        #sns.violinplot(data=df, palette="Set3", bw=.2, cut=1, linewidth=1)
        ## The plot is not over the actual temperatures, the temperatures ara only labels
        #ax.plot(range(len(weights_info[0,:,1])), weights_info[0,:,1], '-o',  label = 'Initial weights')
        #ax.set(
        #    title = 'Weights per state over the entire simulation',
        #    xlabel = 'Sate',
        #    ylabel = 'Weight',
        #)
        #plt.legend()
        #sns.despine(left=True, bottom=True)
        #plt.show()
        ##fig.savefig(os.path.join(dir,'weights_per_state.svg'), bbox_inches="tight")
        #sns.reset_defaults()

    return time, weights_info


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-edr',
                        action='store',
                        dest='edr',
                        type=str)
    parser.add_argument('-log',
                        action='store',
                        dest='log',
                        type=str)
    parser.add_argument('-at',
                        dest='at',
                        nargs='*')
    parser.add_argument('-temp',
                        dest='temp',
                        nargs='*')
    parser.add_argument('-plot',
                        dest='plot',
                        action='store_true')
    args = parser.parse_args()

    #if simulated annealing run
    if args.edr != None:
        #print(args.edr,args.at,args.temp,sep=' ')
        # For get_energy
        annealing_times = [int(s) for s in args.at]
        energy_type = 'Potential'

        # For PandeWeights
        #print(args.temp)
        #print(int(args.temp[0]),int(args.temp[1]+args.temp[2]), int(args.temp[2]))
        temperatures = np.arange(int(args.temp[0]),  int(args.temp[1])+int(args.temp[2]), int(args.temp[2]))
        #print(temperatures)
        energies = get_energy(edr=args.edr, annealing_times=annealing_times, energy_type=energy_type)
        alpha = None
        weight_GROMACS_format = True

        print(energies)
        my_weights = PandeWeights(temperatures=temperatures, energies=energies,alpha=alpha, weight_GROMACS_format=weight_GROMACS_format)
        lambdas = (temperatures-temperatures.min()) / (temperatures.max() - temperatures.min())
        print(*my_weights, *lambdas, sep='\n')
#   if simulating tempering run
#   Getting weight from log
    else:
        time, weights = get_weights_from_log(log=args.log, plot=args.plot)
        fw = weights[-1, :, 1]
        print(*fw)
