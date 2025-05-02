from neuron import h, gui
import neuron
from neuron.units import ms, mV
import time as clock
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from sklearn.metrics import mean_squared_error
import sklearn
import scipy
from scipy.signal import find_peaks
import matplotlib.cm as cm
import math
import pandas as pd
from neuron import coreneuron
import seaborn as sns
import scipy.stats as stats

h.load_file("stdrun.hoc")
h.load_file("import3d.hoc")
h.load_file('nrngui.hoc')
pc = h.ParallelContext()

#### Following functions set up the model
def instantiate_swc(filename):
    ''' 
    Load swc file and instantiate it as cell
    Code source: https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3257
    '''

    # load helper library, included with Neuron
    h.load_file('import3d.hoc')
    # h.load_file(filename)

    # load data
    cell = h.Import3d_SWC_read()
    cell.input(filename)

    # # instantiate
    i3d = h.Import3d_GUI(cell,0)
    i3d.instantiate(None)

def change_Ra(ra=None, electrodeSec=None, electrodeVal=None):
    """Function to change axial resistivity of all model sections"""
    for sec in h.allsec():
        sec.Ra = ra
    if electrodeSec is not None:
        electrodeSec.Ra = electrodeVal

def change_gLeak(gleak= None, erev=None, electrodeSec=None, electrodeVal=None):
    """Function to change leak conductance of all model sections"""
    for sec in h.allsec():
        sec.insert('pas')
        for seg in sec:
            seg.pas.g = gleak
            seg.pas.e = erev
    if electrodeSec is not None:
        for seg in electrodeSec:
            seg.pas.g = electrodeVal
            seg.pas.e = 0

def change_memCap(memcap=1.0, electrodeSec=None, electrodeVal=None):
    """Function to change membrane capacitance all model sections"""
    for sec in h.allsec():
        sec.cm = memcap
    if electrodeSec is not None:
        electrodeSec.cm = electrodeVal

def nsegDiscretization(sectionListToDiscretize):
    """
    This function adjusts the number of segments (nseg) for each section in the given list 
    based on its electrotonic properties. It calculates the electrotonic length constant (λ) 
    and ensures that the segment length does not exceed 10% of λ. If necessary, it increases 
    the number of segments to the next power of 2, rounding up to the nearest odd number 
    to maintain symmetry in discretization. """
    for sec in sectionListToDiscretize:
        secDensityMechDict = sec.psection()['density_mechs']
        secLambda = math.sqrt( ( (1 / secDensityMechDict['pas']['g'][0]) * sec.diam) / (4*sec.Ra) )
        if (sec.L/sec.nseg)/secLambda > 0.1:
            numSeg_log2 = math.log2(sec.L/secLambda / 0.1)
            numSeg = math.ceil(2**numSeg_log2)
            if numSeg % 2 == 0:
                numSeg += 1
            sec.nseg = numSeg

    return

def createElectrode(somaSection, pySectionList):
    """Creation of electrode, and the appended section of the electrode"""
    electrodeSec = h.Section()
    electrodeSec.L = 10
    electrodeSec.diam = 1
    electrodeSec.connect(somaSection, 0)

    pySectionList.append(electrodeSec)

    return pySectionList, electrodeSec

def initializeModel(morph_file, resting_ev=None, sizsection = None, Channels = False, capacitance = None, axial_resistivity = None, gleak_conductance= None, erev = None):
    """Initializes the model with the given parameters:
    resting_ev: membrane potential the cell sits at, which gets changed for simulations when initializing
    sizsection: The spike initiation zone of the model, highlighted in the .swc file with individual nodes labled "dend_12" 
    Channels: False, makes the model passive regardless of what is in this initialize function, (if True, it will give the specified compartments the specified conductances for each channl
    capacitance: Capacitance to give the cell 
    axial_resistivity = axial resitivity of the model
    gleak_conductance= leak conductance of the cell
    erev: reversal potential of the leak channel 
    """
    coreneuron.enable = True
    pc = h.ParallelContext()
    h.cvode.cache_efficient(1)
    cell = instantiate_swc(morph_file)

    allSections_nrn = h.SectionList()
    for sec in h.allsec():
        allSections_nrn.append(sec=sec)
    
    # Create a Python list from this SectionList
    # Select sections from the list by their index

    allSections_py = [sec for sec in allSections_nrn]

    #Define SIZ section
    if sizsection == None:
        sizIndex = 0
    else:
        sizIndex == sizsection
    axonList = h.SectionList()
    tetherList = h.SectionList()
    dendList = h.SectionList()
    preSIZ = h.SectionList()

    for sec in allSections_py:
        if "soma" in sec.name():
            somaSection = sec
        elif "axon" in sec.name():
            axonList.append(sec)
        elif "dend_11" in sec.name():
            tetherList.append(sec)
        elif "dend_12" in sec.name():
            sizSection = sec
        elif "dend_13" in sec.name():
            preSIZ.append(sec)
        else:
            dendList.append(sec)


    allSections_py, electrodeSec = createElectrode(somaSection, allSections_py)

    h.finitialize(resting_ev * mV)
    # h.stdinit()
    if Channels == True:
        #Adding channels throughout the entire morphology 
        for sec in h.allsec():
            sec.insert('na_bac_strege') #The nacbac channel added everywhere 
            for seg in sec:
                seg.gmax_na_bac_strege = 0.00005
                seg.ena = 60
        pct = 1
        for seg in somaSection:
                seg.gmax_na_bac_strege = 0.015
                seg.ena = 60
              
            
        #Adding channels into axon
        h.nat.insert(h.axon)
        h.ks_gunay.insert(h.axon)
        for seg in h.axon:
            seg.natgmax_nat = 0.01875
            seg.ksgmax_ks_gunay = 0.015
            seg.ena = 60
            seg.ek = -80

        #Adding ion channels in SIZ only
        h.nat.insert(sizSection)
        h.nap.insert(sizSection)
        h.ks_gunay.insert(sizSection)

        for seg in sizSection:
            seg.natgmax_nat = 0.03
            seg.napgmax_nap = 0.00018
            seg.ksgmax_ks_gunay = 0.2
            seg.gmax_na_bac_strege = 0.0003125*pct
            seg.ena = 60
            seg.ek = -80
            
    else:
        pass

    change_memCap(memcap=capacitance)
    change_Ra(ra=axial_resistivity)
    change_gLeak(gleak= gleak_conductance, erev=erev)

    nsegDiscretization(allSections_py)
        

    return cell, allSections_py, allSections_nrn, somaSection, sizSection, erev, axonList, tetherList, dendList, electrodeSec

# Function to run the simulation
def run():
    """Run function for utilizing coreneuron"""
    coreneuron.enable = True
    coreneuron.verbose = 0
    coreneuron.model_stats = True
    coreneuron.num_gpus = 1
    pc.set_maxstep(10)
    # h.stdinit()
    pc.barrier()
    pc.psolve(h.tstop)

#####################
#Passive property fitting

def passive_fitting(initial_mV, electrodeSec, somaSection, currents=None, continueRun=200, injDur=None, delay=None, type=None, flynum = None,):
    """
    Parameters:
        initial_mV: resting membrane potential to initialize the model from
        electrodeSec: electrode section inserted into soma
        somaSection: section to insert the electrode into, and stimulate using a current clamp point process
        currents: The current injections you want to provide for the simulations (note this is in nA)
        continueRun: How long the simulation will go for 
        injDur: Duration of the step wise current injection (ms)
        delay: Delay of time before current injections (ms)
        type: Can be "nachbac" or "wt" depending on the cell you are trying to fit
        flynum: The number of which fly data you will be comparing against
    
    Returns:
        simulated_traces: All the traces from the simulations.
        time_points_np: Provides individual time points as a numpy array.
        current_export: Provides a specific current injection simulation used for RMSE calculations, can be changed to be any current injenjecito, but set to -30 pA.
    """
    h.dt = 0.1
    coreneuron.enable = True
    coreneuron.verbose = 0
    coreneuron.model_stats = True
    coreneuron.num_gpus = 1
    h.v_init = initial_mV  # Resting membrane potential

    # Reading in data, change which individual fly you look at here
    if type == "wt":
        data = pd.read_csv(f'Ephys data/Original control wt data/fly_{flynum}_wt_control.csv')
    elif type == "nachbac":
        data = pd.read_csv(f'Ephys data/Original nachbac data/fly_{flynum}_nachbac.csv')

    # Extract time points and membrane potential data
    time_points = data.iloc[:, 0] / 100  
    time_points_np = time_points.values

    # Normalize experimental traces
    hyperpol_current = data.iloc[:, 1:] 
    initial_avg_per_trace = hyperpol_current.iloc[:200, :].mean()
    overall_avg_initial = initial_avg_per_trace.mean()
    normalized_hyperpol_current = hyperpol_current.subtract(initial_avg_per_trace, axis=1).add(overall_avg_initial)
    current_export = normalized_hyperpol_current.iloc[:, 3]  # Adjust what column to compare RMSE

    simulated_traces = []  # Store all traces

    # Run simulations for each current
    if currents is None:
        currents = [0]  # Default to a single current of 0 if none provided

    for current in currents:
        # Set up current clamp
        stimobj = h.IClamp(somaSection(0.5))
        stimobj.delay = delay
        stimobj.dur = injDur
        stimobj.amp = current
        stimobj.i = 0
        ampInjvect = h.Vector().record(stimobj._ref_i)

        vInjVec = h.Vector()
        tInjVec = h.Vector()
        vInjVec.record(somaSection(0.5)._ref_v)
        tInjVec.record(h._ref_t)

        eInjVec = h.Vector()
        eInjVec.record(electrodeSec(0.5)._ref_v)

        h.finitialize(initial_mV * mV)
        h.tstop = continueRun  # ms
        cnargs = coreneuron.nrncore_arg(h.tstop)
        run()

        pc.barrier()
        pc.nrncore_run(cnargs, 1)

        # Retrieve simulation data
        vInj_np = np.array(vInjVec.to_python())
        tInj_np = np.array(tInjVec.to_python())
        eInjVec_np = np.array(eInjVec.to_python())

        # Store each trace
        simulated_traces.append((tInj_np, eInjVec_np, vInj_np, current))

    # Plot experimental vs simulated traces
    fig = plt.figure()
    columns_to_plot = normalized_hyperpol_current.columns[0:5]
    for column in columns_to_plot:
        plt.plot(time_points, normalized_hyperpol_current[column], label=column)
    for tInj_np, eInjVec_np, vInj_np, current in simulated_traces:
        plt.plot(tInj_np, eInjVec_np, label=f'{current*1000} pA', linestyle='-', color='black')
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane potential (mV)')
    plt.legend()
    plt.show()

    pc.psolve(h.tstop)

    return simulated_traces, time_points_np, current_export

def Calculate_RMSE(start_index, end_index, time_data, current_data, simulated_traces, trace_index):
    """
    Compare a selected simulated trace against experimental data.

    Parameters:
        start_index, end_index: Indices defining the range for RMSE calculation.
        time_data, current_data: Experimental time and voltage data.
        simulated_traces: Full list of simulated traces from passive_fitting.
        trace_index: Index of the specific simulated trace to use.

    Returns:
    - RMSE value for the selected trace.
    """
    # Select the desired simulation trace
    tInj_np, eInjVec_np, vInj_np, current = simulated_traces[trace_index]

    # Extract range for experimental data
    time_points_range = time_data[start_index:end_index]
    exp_voltage_range = current_data[start_index:end_index]

    # Interpolate simulation data to match experimental time points
    interpolated_sim_data = np.interp(time_points_range, tInj_np, eInjVec_np)

    # Calculate percent error
    percent_error = np.abs((interpolated_sim_data - exp_voltage_range) / exp_voltage_range) * 100
    average_percent_error = np.mean(percent_error)
    print("Percent Error:", percent_error)
    print("Average Percent Error:", average_percent_error)

    # Calculate RMSE
    rmse = np.sqrt(mean_squared_error(exp_voltage_range, interpolated_sim_data))

    # Plot Experimental vs Selected Simulated Trace
    plt.figure(figsize=(10, 5))
    plt.plot(time_data, current_data, label="Experimental (Full)", color='blue', alpha=0.5)
    plt.plot(tInj_np, eInjVec_np, label=f"Simulated (Full) - {current*1000} pA", linestyle='--', color='gray', alpha=0.5)
    plt.plot(time_points_range, exp_voltage_range, label="Experimental (Selected)", color='blue', linewidth=2)
    plt.plot(time_points_range, interpolated_sim_data, label="Simulated (Selected)", linestyle='--', color='black', linewidth=2)
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane Potential (mV)')
    plt.legend()
    plt.title("Experimental vs Simulated Voltage Traces")
    plt.show()

    print("RMSE:", rmse)
    return rmse

#####################
#Helper functions
def normalize_traces(type = None, flynum = None, plots = False, traces = False, current_start = None, current_end = None):
    """ There are 2 possible types wt, nachbac
    there are 10 flies for wt
    There are 14 flies for nachbac
    Function base line subtracts each of traces to the average resting membrane potential using the initial 200 ms of the recording prior to injection any current

    Parameters:
        type: Can be "nachbac" or "wt" depending on the cell you are trying to fit
        flynum: The number of which fly data you will be comparing against
        plots: If False will just return the data (hyperpolarizing current injection and time points)
        traces: If false then it will use current_start and current_end that are provided, if True defaults to -30 to 60
        current_start: The beginning current injection to normalize from (pA)
        current_end: The ending current injection to normalize (pA)
    Returns:
        Time_points: a data column with indidivual time points, and normalized_current the current injections that have been baseline subtracted
    """
    flynum = str(flynum) 
    if type == "wt":
        data = pd.read_csv(f'Ephys data/Original control wt data/fly_{flynum}_wt_control.csv')
    elif type == "nachbac":
        data = pd.read_csv(f'Ephys data/Original nachbac data/fly_{flynum}_nachbac.csv')
    if traces == False:
        current_start = current_start
        current_end = current_end
    elif traces == True:
        current_start = -30
        current_end = 60
    current_value = current_start

    # Extract time points and hyperpolarization current data
    time_points = data.iloc[:, 0] / 100  # First column, divided by 100 (convert to seconds)
    current_injections = data.iloc[:, 1:]  # Remaining columns are membrane potentials

    # Calculate the average membrane potential of the initial points
    hyperpol_current = data.iloc[:, 1:]
    
    initial_avg_per_trace = hyperpol_current.iloc[:200, :].mean()

    # Calculate the overall average of these initial 200 ms means
    overall_avg_initial = initial_avg_per_trace.mean()

    # Normalize each trace: Subtract its own initial 200 ms mean and add the overall average
    normalized_current = hyperpol_current.subtract(initial_avg_per_trace, axis=1).add(overall_avg_initial)

    # Plotting the original and normalized traces side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    start_column_index = (current_start + 40) // 10
    end_column_index = (current_end + 40) // 10
    current_start = current_start
    for column in current_injections.columns[start_column_index:end_column_index + 1]:
        # Plot the trace with the current value as the label
        ax1.plot(time_points, current_injections[column], label=f'{current_value} pA')
        
        # Increment the current value for the label
        current_value += 10

    # Set labels and title
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Membrane Potential (mV)')
    ax1.set_title(f'Original Current Traces from {type} fly {flynum}')
    # ax1.legend(loc='upper right', bbox_to_anchor=(1.05, 1), fontsize='small', frameon=False)

# Plot normalized traces
    start_column_index = (current_start + 40) // 10
    end_column_index = (current_end + 40) // 10

    # Plot normalized traces only within the specified current range
    current_value = current_start
    for column in normalized_current.columns[start_column_index:end_column_index + 1]:
        # Plot the trace with the current value as the label
        ax2.plot(time_points, normalized_current[column], label=f'{current_value} pA')
        
        # Increment the current value for the label
        current_value += 10

    # Set labels and title
    ax2.set_xlabel('Time (ms)')
    ax2.set_ylabel('Membrane Potential (mV)')
    ax2.set_title(f'Current Traces from {type} fly {flynum}')
    # ax2.legend(loc='upper right', bbox_to_anchor=(1.05, 1), fontsize='small', frameon=False)

    handles, labels = ax1.get_legend_handles_labels()  # Get handles and labels from one of the subplots
    fig.legend(handles, labels, loc='lower center', fontsize='medium', frameon=False, ncol=len(labels) // 2)  # Arrange in two rows

    plt.tight_layout(rect=[0, 0.1, 1, 1])
    if plots == True:
        plt.show()
    elif plots == False:
        return time_points, normalized_current

def plot_iv_curve(simulated_traces, injDur, delay, inject_time):
    """
    Generates an I-V curve by plotting the peak current against the membrane potential.
    
    Parameters:
        simulated_traces: List of dictionaries containing simulation results.
                          Each dictionary has 'time', 'current', 'voltage', and 'step'.
        injDur: Duration of the injection (ms).
        delay: Delay before injection (ms).
    """
    peak_currents = []
    steady_state_currents = []
    membrane_potentials = []
    
    for trace in simulated_traces:
        time = trace['time']
        current = trace['current']
        step = trace['step']  # Membrane potential
        
        # Identify indices for the injection period
        start_idx = np.where(time >= +inject_time)[0][0]  # Start of injection
        start_idx2 = np.where(time >= inject_time+0.2)[0][0]  # Start of injection
        end_idx = np.where(time >= (inject_time + injDur))[0][0]  # End of injection
        
        # Calculate steady-state current: Mean of the last 2 ms of the injection
        steady_state_idx = np.where(time >= (inject_time + injDur - 20))[0]
        steady_state_current = np.mean(current[steady_state_idx])
        
        # Calculate peak current
        peak_current = np.min(current[start_idx2:end_idx]) - steady_state_current
        
        # Append results
        steady_state_currents.append(steady_state_current)
        peak_currents.append(peak_current*100)
        membrane_potentials.append(step)
    
    # Plot the I-V curve
    plt.figure(figsize=(8, 6))
    plt.plot(membrane_potentials, peak_currents, 'o-', label="Peak Current")
    # plt.plot(membrane_potentials, steady_state_currents, 's--', label="Steady-State Current")
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
    plt.title("I-V Curve")
    plt.xlabel("Membrane Potential (mV)", fontsize =12)
    plt.ylim(-100, 10)
    plt.ylabel("Current (pA)", fontsize =12)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_spike_count_vs_current(csv_file):
    """
    Plots Spike Count vs. Current for each unique conductance, sorted from lowest to highest.

    Parameters:
        csv_file (str): Path to the CSV file.
    """
    # Read the CSV file
    data = pd.read_csv(csv_file)

    # Ensure the necessary columns exist
    required_columns = {'Current (nA)', 'Spike Count', 'Conductance of nachbac'}
    if not required_columns.issubset(data.columns):
        raise ValueError(f"CSV file must contain columns: {required_columns}")

    # Get unique conductance values, sorted in ascending order
    unique_conductances = sorted(data['Conductance of nachbac'].unique())

    # Set seaborn style
    sns.set_theme(style="whitegrid")

    # Generate distinct colors for each conductance value
    colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_conductances)))

    # Create the plot
    plt.figure(figsize=(10, 6))

    for i, conductance in enumerate(unique_conductances):
        # Filter data for the current conductance value
        subset = data[data['Conductance of nachbac'] == conductance]

        # Sort the subset by 'Current (nA)' for consistent plotting
        subset = subset.sort_values(by='Current (nA)')

        # Plot with unique color
        plt.plot(
            subset['Current (nA)'], 
            subset['Spike Count'], 
            marker='o', 
            color=colors[i], 
            label=f'Conductance = {conductance}'
        )

    # Add labels, title, and legend
    plt.xlabel('Current (nA)')
    plt.ylabel('Spike Count')
    plt.title('Spike Count vs Current for Different Conductances')
    plt.legend(title='Conductance of nachbac', loc='best')
    plt.tight_layout()
    plt.show()

def plot_spike_counts_vs_currents(csv1=None, csv2=None):
    """Takes in 2 csv files, WT and nachbac data and generates various plots to visualize spike counts vs current across simulations"""
    # Load data
    wt = pd.read_csv(csv1)
    data = wt.copy()  # Start with WT data
    if csv2 is not None:
        nachbac = pd.read_csv(csv2)
        data = pd.concat([wt, nachbac], ignore_index=True)

    plt.figure(figsize=(12, 6))

    # Set to track which Flynum and Sim_type combinations have been added to the legend
    legend_labels = set()

    # Loop through each unique combination of Flynum and Sim_type and plot
    for (flynum, sim_type) in data.groupby(['Flynum', 'Sim_type', "Para decreased"]):
        fly_data = sim_type
        # Convert Flynum and Sim_type to strings for the legend
        label = f"Flynum {str(flynum)}, Sim_type {str(sim_type['Sim_type'].iloc[0])}"

        # Only add the label if it hasn't been added before
        if label not in legend_labels:
            legend_labels.add(label)
            plt.plot(fly_data['Current (nA)'], fly_data['Spike Count'], marker='o', label=label)
        else:
            plt.plot(fly_data['Current (nA)'], fly_data['Spike Count'], marker='o')

    plt.title("Spike Count vs Current for Each Flynum and Sim_type")
    plt.xlabel("Current (nA)")
    plt.ylabel("Spike Count")
    plt.legend(title="Flynum and Sim_type", loc="upper left")
    plt.tight_layout()
    plt.show()

     # Create a figure for Para decreased = False
    plt.figure(figsize=(12, 6))
    plt.title("WT vs NaChBac models")
    legend_labels = set()
    for (flynum, sim_type) in data[data['Para decreased'] == False].groupby(['Flynum', 'Sim_type', 'Para decreased']):
        fly_data = sim_type
        label = f"Flynum {str(flynum)}, Sim_type {str(sim_type['Sim_type'].iloc[0])}"

        if label not in legend_labels:
            legend_labels.add(label)
            plt.plot(fly_data['Current (nA)'], fly_data['Spike Count'], marker='o', label=label)
        else:
            plt.plot(fly_data['Current (nA)'], fly_data['Spike Count'], marker='o')

    plt.xlabel("Current (nA)")
    plt.ylabel("Spike Count")
    plt.legend(title="Flynum and Sim_type", loc="upper left")
    plt.tight_layout()
    plt.show()

    # Create a figure for Para decreased = True
    plt.figure(figsize=(12, 6))
    plt.title("NaChBac models (Para not decreased)")
    legend_labels = set()
    for (flynum, sim_type) in data[data['Para decreased'] == True].groupby(['Flynum', 'Sim_type', 'Para decreased']):
        fly_data = sim_type
        label = f"Flynum {str(flynum)}, Sim_type {str(sim_type['Sim_type'].iloc[0])}"

        if label not in legend_labels:
            legend_labels.add(label)
            plt.plot(fly_data['Current (nA)'], fly_data['Spike Count'], marker='o', label=label)
        else:
            plt.plot(fly_data['Current (nA)'], fly_data['Spike Count'], marker='o')

    plt.xlabel("Current (nA)")
    plt.ylabel("Spike Count")
    plt.legend(title="Flynum and Sim_type", loc="upper left")
    plt.tight_layout()
    plt.show()

    para_false_data = data[data['Para decreased'] == False]

    avg_data = para_false_data.groupby(['Sim_type', 'Current (nA)']).agg(
        avg_spike_count=('Spike Count', 'mean')
    ).reset_index()

    # Create the plot
    plt.figure(figsize=(12, 6))
    plt.title("WT vs NaChBac models")

    # Set to track which Sim_type combinations have been added to the legend
    legend_labels = set()

    # Loop through each unique Sim_type and plot
    for sim_type in avg_data['Sim_type'].unique():
        sim_data = avg_data[avg_data['Sim_type'] == sim_type]
        label = f"Sim_type {str(sim_type)}"

        if label not in legend_labels:
            legend_labels.add(label)
            plt.plot(sim_data['Current (nA)'], sim_data['avg_spike_count'], marker='o', label=label)
        else:
            plt.plot(sim_data['Current (nA)'], sim_data['avg_spike_count'], marker='o')

    plt.xlabel("Current (nA)")
    plt.ylabel("Average Spike Count")
    plt.legend(title="Sim_type", loc="upper left")
    plt.tight_layout()
    plt.show()

    para_True_data = data[data['Para decreased'] == True]

    avg_data = para_True_data.groupby(['Sim_type', 'Current (nA)']).agg(
        avg_spike_count=('Spike Count', 'mean')
    ).reset_index()

    # Create the plot
    plt.figure(figsize=(12, 6))
    plt.title("NaChBac models (para not decreased)")

    # Set to track which Sim_type combinations have been added to the legend
    legend_labels = set()

    # Loop through each unique Sim_type and plot
    for sim_type in avg_data['Sim_type'].unique():
        sim_data = avg_data[avg_data['Sim_type'] == sim_type]
        label = f"Sim_type {str(sim_type)}"

        if label not in legend_labels:
            legend_labels.add(label)
            plt.plot(sim_data['Current (nA)'], sim_data['avg_spike_count'], marker='o', label=label)
        else:
            plt.plot(sim_data['Current (nA)'], sim_data['avg_spike_count'], marker='o')

    plt.xlabel("Current (nA)")
    plt.ylabel("Average Spike Count")
    plt.legend(title="Sim_type", loc="upper left")
    plt.tight_layout()
    plt.show()

    wt_data = data[data['Sim_type'] == 'WT']
    nachbac_data = data[data['Sim_type'] != 'WT']

    # Filter WT data for Para decreased = False and True
    wt_para_false = wt_data[wt_data['Para decreased'] == False]
    wt_para_true = wt_data[wt_data['Para decreased'] == True]

    # Filter nachbac data for Para decreased = False and True
    nachbac_para_false = nachbac_data[nachbac_data['Para decreased'] == False]
    nachbac_para_true = nachbac_data[nachbac_data['Para decreased'] == True]

    # Group by Current (nA) and calculate average Spike Count and standard deviation
    wt_avg_data_false = wt_para_false.groupby('Current (nA)').agg(
        avg_spike_count=('Spike Count', 'mean'),
        std_spike_count=('Spike Count', 'std')
    ).reset_index()

    wt_avg_data_true = wt_para_true.groupby('Current (nA)').agg(
        avg_spike_count=('Spike Count', 'mean'),
        std_spike_count=('Spike Count', 'std')
    ).reset_index()

    nachbac_avg_data_false = nachbac_para_false.groupby('Current (nA)').agg(
        avg_spike_count=('Spike Count', 'mean'),
        std_spike_count=('Spike Count', 'std')
    ).reset_index()

    nachbac_avg_data_true = nachbac_para_true.groupby('Current (nA)').agg(
        avg_spike_count=('Spike Count', 'mean'),
        std_spike_count=('Spike Count', 'std')
    ).reset_index()

    # Define a helper function to ensure error bars don't extend below zero
    def plot_error_bars(x, y, yerr, **kwargs):
        lower_error = np.minimum(yerr, y)  # Limit lower error bar to the value of y
        upper_error = yerr
        plt.errorbar(x, y, yerr=[lower_error, upper_error], **kwargs)

    # Create the plot
    plt.figure(figsize=(12, 6))
    plt.title("Average Spike Count vs Current (WT and nachbac, Para decreased = False and True)")

    # Plot for WT, Para decreased = False
    plt.plot(wt_avg_data_false['Current (nA)'], wt_avg_data_false['avg_spike_count'], 
            'o-', label="WT")  # Solid line and points
    plot_error_bars(wt_avg_data_false['Current (nA)'], wt_avg_data_false['avg_spike_count'], 
                    wt_avg_data_false['std_spike_count'], fmt='none', ecolor='blue', alpha=0.5, capsize=5)

    # Plot for nachbac, Para decreased = False
    plt.plot(nachbac_avg_data_false['Current (nA)'], nachbac_avg_data_false['avg_spike_count'], 
            'o-', label="NaChBac, Para decreased")  # Solid line and points
    plot_error_bars(nachbac_avg_data_false['Current (nA)'], nachbac_avg_data_false['avg_spike_count'], 
                    nachbac_avg_data_false['std_spike_count'], fmt='none', ecolor='orange', alpha=0.5, capsize=5)

    # Plot for nachbac, Para decreased = True
    plt.plot(nachbac_avg_data_true['Current (nA)'], nachbac_avg_data_true['avg_spike_count'], 
            'o-', label="NaChBac, Para not decreased")  # Solid line and points
    plot_error_bars(nachbac_avg_data_true['Current (nA)'], nachbac_avg_data_true['avg_spike_count'], 
                    nachbac_avg_data_true['std_spike_count'], fmt='none', ecolor='green', alpha=0.5, capsize=5)

    # Add labels and legend
    plt.xlabel("Current (nA)")
    plt.ylabel("Average Spike Count")
    plt.legend(title="Sim_type and Para decreased", loc="upper left")
    plt.tight_layout()
    plt.show()

def NaChBac_tau_plotting():
    strege_taum = pd.read_csv('Strege et al tau.csv')
    strege_tauh = pd.read_csv('Strege et al tauh.csv')

    def tauh_strege(v):
        return 36.11944425 / ( 1 + np.exp((-24.58946004-v)/-5.15842944)) + 165.144429082
      
    def taum_strege(v):
        return -121675.59 / (1 + np.exp(-0.06 * (v - -183.26))) + 121680.09

        
    v_values = np.linspace(-50, 30, 400)

    m_values_strege = taum_strege(v_values)
    h_values_strege = tauh_strege(v_values)


    plt.figure(figsize=(10, 6))

    # Plot m vs v
    plt.subplot(2, 1, 1)
    plt.plot(v_values, m_values_strege, label='taum(v) strege', color='green')
    plt.errorbar(strege_taum['Membrane potential'], strege_taum['taum'], yerr=strege_taum["SEM"], fmt='o', color='black', capsize=3, label='taum(v) values strege')
    plt.xlabel('v (mV)')
    plt.ylabel('Tau (activation)')
    plt.title('Time constant over voltage')
    plt.grid(True)
    plt.legend()

  

    plt.subplot(2, 1, 2)
    plt.plot(v_values, h_values_strege, label='tauh(v) strege', color='green')
    plt.errorbar(strege_tauh['Membrane potential'], strege_tauh['tauh'], yerr=strege_tauh["SEM"], fmt='o', color='black', capsize=3, label='tauh(v) values strege')
    plt.xlabel('v (mV)')
    plt.ylabel('Tau (inactivation)')
    plt.grid(True)
    plt.legend()

    plt.show()


#####################

#Simulations

def Varying_nachbac_soma_current_injections(initial_mV, electrodeSec, somaSection, currents=None, continueRun=None, injDur=None, delay = None, type = None):
    """"Simulation of current injections into soma of model, keeps track of the current, nachbac cond. and the number of spikes for each injection
    Parameters:
        initial_mV: resting membrane potential to initialize the model from
        electrodeSec: electrode section inserted into soma
        somaSection: section to insert the electrode into, and stimulate using a current clamp point process
        currents: The current injections you want to provide for the simulations (note this is in nA)
        continueRun: How long the simulation will go for 
        injDur: Duration of the step wise current injection (ms)
        delay: Delay of time before current injections (ms)
        type: Can be "nachbac" or "wt" depending on the cell you are trying to fit.

    Returns:
        -saves a csv file with the data Current (nA), Spike Count Firing Frequency (Hz), Conductance of nachbac
    }


    """
    h.dt = 0.1
    print(h.dt)
    coreneuron.enable = True
    coreneuron.verbose = 0
    coreneuron.model_stats = True
    coreneuron.num_gpus = 1
    h.v_init = initial_mV

    simulated_traces = []

    # Ask for percentage of Nat/Nap
    if type == "nachbac":
        data_column = float(input("Conductance of Nacbac: "))


    # Set up lists to store analysis results
    spike_counts = []
    firing_frequencies = []
    currents_used = []

    # Set the minimum height for a peak relative to the surrounding depolarization
    relative_spike_threshold = 2  # in mV, adjustable based on expected spike amplitude

    # Default currents if not provided
    if currents is None:
        currents = [0]  # Default to a single current of 0 if none provided

    # Run simulations for each current
    for current in currents:
        # Set up current clamp
        stimobj = h.IClamp(somaSection(0.5))
        stimobj.delay = delay
        stimobj.dur = injDur
        stimobj.amp = current
        stimobj.i = 0
        ampInjvect = h.Vector().record(stimobj._ref_i)

        vInjVec = h.Vector()
        tInjVec = h.Vector()
        vInjVec.record(somaSection(0.5)._ref_v)
        tInjVec.record(h._ref_t)

        eInjVec = h.Vector()
        eInjVec.record(electrodeSec(0.5)._ref_v)

        # Initialize and run the simulation
        h.finitialize(initial_mV * mV)
        h.tstop = continueRun  # ms
        cnargs = coreneuron.nrncore_arg(h.tstop)
        run()

        # Synchronize across processors
        pc.barrier()
        pc.nrncore_run(cnargs, 1)

        # Retrieve simulation data
        aInj = ampInjvect.to_python()
        vInj = vInjVec.to_python()
        tInj = tInjVec.to_python()
        vInj_np = np.array(vInj)
        tInj_np = np.array(tInj)

        # Store simulation data
        simulated_traces.append((tInj_np, vInj_np, current))

        # Detect spikes based on relative threshold
        peaks, _ = find_peaks(vInj_np, height=(np.mean(vInj_np) + relative_spike_threshold))
        spike_times = [tInj_np[p] for p in peaks if stimobj.delay <= tInj_np[p] <= stimobj.delay + stimobj.dur and tInj_np[p] >= 630]
        spike_count = len(spike_times)
        duration_in_seconds = injDur / 1000.0  # Convert ms to seconds
        firing_frequency = spike_count / duration_in_seconds if duration_in_seconds > 0 else 0

        # Store the results
        spike_counts.append(spike_count)
        firing_frequencies.append(firing_frequency)
        currents_used.append(current)

    plt.figure()
    for tInj_np, vInj_np, current in simulated_traces:
        plt.plot(tInj_np, vInj_np, label=f'Current: {current} nA')

    # Labeling the plot
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane Potential (mV)')
    plt.title('Current injection in Somas')
    plt.legend()
    plt.grid()
    plt.show()

    # Plot spike count vs. current
    plt.figure()
    plt.plot(currents_used, spike_counts, marker='o')
    plt.xlabel('Injected Current (nA)')
    plt.ylabel('Number of Spikes')
    plt.title('Spike Count vs. Injected Current')
    plt.grid()
    plt.show()

    # Plot firing frequency vs. current
    plt.figure()
    plt.plot(currents_used, firing_frequencies, marker='o')
    plt.xlabel('Injected Current (nA)')
    plt.ylabel('Firing Frequency (Hz)')
    plt.title('Firing Frequency vs. Injected Current')
    plt.grid()
    plt.show()

    # Save data to CSV
    if type == "nachbac":
         data = {
        'Current (nA)': currents_used,
        'Spike Count': spike_counts,
        'Firing Frequency (Hz)': firing_frequencies,
        'Conductance of nachbac': [data_column] * len(currents_used)
    }
    else:
        data = {
            'Current (nA)': currents_used,
            'Spike Count': spike_counts,
            'Firing Frequency (Hz)': firing_frequencies,
            'Percentage Nat/Nap': [data_column] * len(currents_used)
        }
    df = pd.DataFrame(data)
   # Save data to CSV, appending if the file exists
    if type == "nachbac":
        csv_filename = 'WT_simulation_with_varying_nachbac_cond.csv'
    else:     
        csv_filename = 'WT_simulation_results2.csv'
    if os.path.exists(csv_filename):
        # If file exists, load the existing data and append the new data
        existing_df = pd.read_csv(csv_filename)
        combined_df = pd.concat([existing_df, df], ignore_index=True)
        combined_df.to_csv(csv_filename, index=False)
        print(f"New simulation results appended to {csv_filename}")
    else:
        # If file does not exist, create a new file
        df.to_csv(csv_filename, index=False)
        print(f"Simulation results saved to {csv_filename}")

def Soma_current_injections(initial_mV, electrodeSec, somaSection, currents=None, continueRun=None, injDur=None, delay = None, type = None, manipulation = False, flynum = None):
    """"Simulation is a sensitivity analysis, simple current injection simulation with given current in pA across conditions WT, and  nachbac
    This is to compare the firing frequencies across model types for giving current injections.
    
    Parameters:
        initial_mV: resting membrane potential to initialize the model from
        electrodeSec: electrode section inserted into soma
        somaSection: section to insert the electrode into, and stimulate using a current clamp point process
        currents: The current injections you want to provide for the simulations (note this is in nA)
        continueRun: How long the simulation will go for 
        injDur: Duration of the step wise current injection (ms)
        delay: Delay of time before current injections (ms)
        type: Can be "nachbac" or "wt" depending on the cell you are trying to fit.
        manipulation: If true then para was NOT decreased in the simulation, if false then para was decreased
        flynum: Which fly the simulation is associated with (which passive property set is used)

    Returns:
    - a csv with currents (nA), spike Counts, firing frequency (Hz) sim_type, whether Para was decreased or not, and the fly the simulation is associated with
    }
    
    """
    
    h.dt = 0.1
    print(h.dt)
    coreneuron.enable = True
    coreneuron.verbose = 0
    coreneuron.model_stats = True
    coreneuron.num_gpus = 1
    h.v_init = initial_mV

    simulated_traces = []

    # Ask for the simulation type
    if type == "nachbac" and manipulation == True:
        data_column = input("What NachBac sim is this?: ")
        para_decreased = True

    if type == "nachbac" and manipulation == False:
        data_column = input("What NachBac sim is this?: ")
        para_decreased = False

    if type == "wt":
        data_column = "WT"
        para_decreased = False

    # Set up lists to store analysis results
    spike_counts = []
    firing_frequencies = []
    currents_used = []

    # Set the minimum height for a peak relative to the surrounding depolarization
    relative_spike_threshold = 2  # in mV, adjustable based on expected spike amplitude

    # Default currents if not provided
    if currents is None:
        currents = [0]  # Default to a single current of 0 if none provided

    # Run simulations for each current
    for current in currents:
        # Set up current clamp
        stimobj = h.IClamp(somaSection(0.5))
        stimobj.delay = delay
        stimobj.dur = injDur
        stimobj.amp = current
        stimobj.i = 0
        ampInjvect = h.Vector().record(stimobj._ref_i)

        vInjVec = h.Vector()
        tInjVec = h.Vector()
        vInjVec.record(somaSection(0.5)._ref_v)
        tInjVec.record(h._ref_t)

        eInjVec = h.Vector()
        eInjVec.record(electrodeSec(0.5)._ref_v)

        # Initialize and run the simulation
        h.finitialize(initial_mV * mV)
        h.tstop = continueRun  # ms
        cnargs = coreneuron.nrncore_arg(h.tstop)
        run()

        # Synchronize across processors
        pc.barrier()
        pc.nrncore_run(cnargs, 1)

        # Retrieve simulation data
        aInj = ampInjvect.to_python()
        vInj = vInjVec.to_python()
        tInj = tInjVec.to_python()
        vInj_np = np.array(vInj)
        tInj_np = np.array(tInj)

        # Store simulation data
        simulated_traces.append((tInj_np, vInj_np, current))

        # Detect spikes based on relative threshold
        peaks, _ = find_peaks(vInj_np, height=(np.mean(vInj_np) + relative_spike_threshold))
        if type == "nachbac":
            spike_times = [tInj_np[p] for p in peaks if stimobj.delay <= tInj_np[p] <= stimobj.delay + stimobj.dur and tInj_np[p] >= 630]
        else:
            spike_times = [tInj_np[p] for p in peaks if stimobj.delay <= tInj_np[p] <= stimobj.delay + stimobj.dur]
        spike_count = len(spike_times)
        duration_in_seconds = injDur / 1000.0  # Convert ms to seconds
        firing_frequency = spike_count / duration_in_seconds if duration_in_seconds > 0 else 0

        # Store the results
        spike_counts.append(spike_count)
        firing_frequencies.append(firing_frequency)
        currents_used.append(current)

    plt.figure()
    for tInj_np, vInj_np, current in simulated_traces:
    # Plot the voltage trace
        plt.plot(tInj_np, vInj_np, label=f'Current: {current} nA')

        # Detect spikes using relative threshold
        peaks, _ = find_peaks(vInj_np, height=(np.mean(vInj_np) + relative_spike_threshold))
        if type == "nachbac":
            spike_times = [tInj_np[p] for p in peaks if stimobj.delay <= tInj_np[p] <= stimobj.delay + stimobj.dur and tInj_np[p] >= 630]
        else:
            spike_times = [tInj_np[p] for p in peaks if stimobj.delay <= tInj_np[p] <= stimobj.delay + stimobj.dur]

        # Plot a single marker for each detected spike
        for spike_time in spike_times:
            # Use a tolerance to find the nearest index
            spike_index = np.argmin(np.abs(tInj_np - spike_time))  # Find the closest index
            plt.scatter(tInj_np[spike_index], vInj_np[spike_index], color='red', marker='o' )

    # Labeling the plot
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane Potential (mV)')
    plt.title('Current injection in Somas')
    plt.legend()
    plt.grid()
    plt.show()

    # Plot spike count vs. current
    plt.figure()
    plt.plot(currents_used, spike_counts, marker='o')
    plt.xlabel('Injected Current (nA)')
    plt.ylabel('Number of Spikes')
    plt.title('Spike Count vs. Injected Current')
    plt.grid()
    plt.show()

    # Plot firing frequency vs. current
    plt.figure()
    plt.plot(currents_used, firing_frequencies, marker='o')
    plt.xlabel('Injected Current (nA)')
    plt.ylabel('Firing Frequency (Hz)')
    plt.title('Firing Frequency vs. Injected Current')
    plt.grid()
    plt.show()

    # Save data to CSV
    if type == "nachbac" and manipulation == True :
         data = {
        'Current (nA)': currents_used,
        'Spike Count': spike_counts,
        'Firing Frequency (Hz)': firing_frequencies,
        'Sim_type': [data_column] * len(currents_used),
        'Para decreased': [para_decreased] * len(currents_used),
        'Flynum': flynum
    }
         
    elif type == "nachbac" and manipulation == False :
         data = {
        'Current (nA)': currents_used,
        'Spike Count': spike_counts,
        'Firing Frequency (Hz)': firing_frequencies,
        'Sim_type': [data_column] * len(currents_used),
        'Para decreased': [para_decreased] * len(currents_used),
        'Flynum': flynum
    }
    elif type == "wt":
         data = {
        'Current (nA)': currents_used,
        'Spike Count': spike_counts,
        'Firing Frequency (Hz)': firing_frequencies,
        'Sim_type': [data_column] * len(currents_used),
        'Para decreased': [para_decreased] * len(currents_used),
        'Flynum': flynum
    }

    df = pd.DataFrame(data)
   # Save data to CSV, appending if the file exists
    if type == "nachbac":
        csv_filename = 'simulation_results_nachbac_final.csv'
    elif type == 'wt':     
        csv_filename = 'wt_simulation_results.csv'
    if os.path.exists(csv_filename):
        # If file exists, load the existing data and append the new data
        existing_df = pd.read_csv(csv_filename)
        combined_df = pd.concat([existing_df, df], ignore_index=True)
        combined_df.to_csv(csv_filename, index=False)
        print(f"New simulation results appended to {csv_filename}")
    else:
        # If file does not exist, create a new file
        df.to_csv(csv_filename, index=False)
        print(f"Simulation results saved to {csv_filename}")

def runSimAndNormalize(initial_mV, electrodeSec, somaSection, currents=None, continueRun=200, injDur=None, delay=None, type=None, flynum = None, current_start = None, current_end = None):
    """
    Parameters: 
        initial_mV: resting membrane potential to initialize the model from
        electrodeSec: electrode section inserted into soma
        somaSection: section to insert the electrode into, and stimulate using a current clamp point process
        currents: The current injections you want to provide for the simulations (note this is in nA)
        continueRun: How long the simulation will go for 
        injDur: Duration of the step wise current injection (ms)
        delay: Delay of time before current injections (ms)
        type: Can be "nachbac" or "wt" depending on the cell you are trying to fit.
        flynum: Which fly the simulation is associated with (which passive property set is used)
        current_start: The beginning current injection to normalize from (pA)
        current_end: The ending current injection to normalize (pA)

    """
    
    
    start = clock.time()

    # CoreNEURON settings
    h.dt = 0.1
    print(h.dt)
    coreneuron.enable = True
    coreneuron.verbose = 0
    coreneuron.model_stats = True
    coreneuron.num_gpus = 1
    h.v_init = initial_mV  # Resting membrane potential

    time_points, normalized_hyperpol_current =normalize_traces(type = type, flynum = flynum, plots = False, traces = True, current_start = current_start, current_end = current_end)
    
    # List to store simulation data for each current
    simulated_traces = []

    # Run simulations for each current
    if currents is None:
        currents = [0]  # Default to a single current of 0 if none provided

    for current in currents:
        # Set up current clamp
        stimobj = h.IClamp(somaSection(0.5))
        stimobj.delay = delay
        stimobj.dur = injDur
        stimobj.amp = current
        stimobj.i = 0
        ampInjvect = h.Vector().record(stimobj._ref_i)

        vInjVec = h.Vector()
        tInjVec = h.Vector()
        vInjVec.record(somaSection(0.5)._ref_v)
        tInjVec.record(h._ref_t)

        eInjVec = h.Vector()
        eInjVec.record(electrodeSec(0.5)._ref_v)

        # h.stdinit()
        h.finitialize(initial_mV * mV)
        # Run the simulation
        h.tstop = continueRun  # ms
        cnargs = coreneuron.nrncore_arg(h.tstop)
        run()

        # Arguments for CoreNEURON run
        pc.barrier()
        pc.nrncore_run(cnargs, 1)


        # Retrieve simulation data
        aInj = ampInjvect.to_python()
        vInj = vInjVec.to_python()
        tInj = tInjVec.to_python()
        vInj_np = np.array(vInj)
        tInj_np = np.array(tInj)
        eInjVec_np = np.array(eInjVec)

        # Store simulation data
        simulated_traces.append((tInj_np, eInjVec_np, vInj_np, current))

    currents = currents
    colors = cm.viridis(np.linspace(0, 1, len(currents)))  # Generate colors from a colormap

    # Create a color dictionary to map each current to a unique color
    color_dict = {currents[i]: colors[i] for i in range(len(currents))}

    # Initialize a new figure
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

    # Plot normalized experimental traces
    num_columns = len(normalized_hyperpol_current.columns)
    num_currents = len(currents)

    # Check if the number of columns and currents match
    if num_columns > num_currents:
        print("Warning: More experimental traces than currents. Some traces might not have corresponding currents.")
    elif num_columns < num_currents:
        print("Warning: More currents than experimental traces. Some currents will not be represented.")


    start_column_index = (current_start + 40) // 10
    end_column_index = (current_end + 40) // 10
    current_start = current_start
    current_value = current_start
    for column in normalized_hyperpol_current.columns[start_column_index:end_column_index + 1]:
        # Plot the trace with the current value as the label
        ax1.plot(time_points, normalized_hyperpol_current[column], label=f'{current_value} pA')
        ax1.legend(loc='upper right', fontsize='small', frameon=False)
        
        # Increment the current value for the label
        current_value += 10

    # Set labels and title
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Membrane Potential (mV)')
    ax1.set_title(f'Original Current Traces from {type} fly {flynum}')
    ax1.set_ylim(-130, 60)

    # Plot simulated traces with the same color as their corresponding experimental traces
    for tInj_np, eInjVec_np, vInj_np, current in simulated_traces:
        if current in color_dict:  # Use color from the color dictionary
            ax2.plot(tInj_np, eInjVec_np, label=f'{current*1000} pA', linestyle='-')
        #plt.plot(tInj_np, vInj_np, label=f'Simulated (Soma) - {current} pA', linestyle='-')
            ax2.set_xlabel('Time (ms)')
            ax2.set_ylabel('Membrane Potential (mV)')
            ax2.set_title('Simulated Traces')
            ax2.legend(loc='upper right', fontsize='small', frameon=False)
            ax2.set_ylim(-130, 50)

    for tInj_np, eInjVec_np, vInj_np, current in simulated_traces:
        if current in color_dict:  # Use color from the color dictionary
            ax3.plot(tInj_np, eInjVec_np, label=f'{current*1000} pA', linestyle='-')

    for column in normalized_hyperpol_current.columns[start_column_index:end_column_index + 1]:
        # Plot the trace with the current value as the label
        ax3.plot(time_points, normalized_hyperpol_current[column], label=f'{current_value} pA', color="black")
        
        # Increment the current value for the label
        current_value += 10

    # Set labels and title
    ax3.set_xlabel('Time (ms)')
    ax3.set_ylabel('Membrane Potential (mV)')
    ax3.set_title('Overlayed Traces')
    plt.show()

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    selected_values = [-30, -20, -10, 0, 10, 20, 30, 40, 50, 60]

    # Find the indices of these columns based on the step of +10 starting from -40
    start_value = -40
    step = 10

    # Compute the indices
    selected_indices = [(value - start_value) // step for value in selected_values]

    # Subset the data frame based on these indices
    selected_columns = normalized_hyperpol_current.iloc[:, selected_indices]

    current_value = selected_values[0]  # Initialize to the first value in selected_values

    for column in selected_columns.columns:
        # Plot the trace with the current value as the label
        ax1.plot(time_points, normalized_hyperpol_current[column], label=f'{current_value} pA')
        ax1.legend(loc='upper right', fontsize='small', frameon=False)
        ax3.plot(time_points, normalized_hyperpol_current[column], label=f'{current_value} pA', color="black")
        
        # Increment the current value for the label
        current_value += 10

    # Set labels and title
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Membrane Potential (mV)')
    ax1.set_title(f'Original Current Traces from {type} fly {flynum}')
    ax1.set_ylim(-130, 50)

    for tInj_np, eInjVec_np, vInj_np, current in simulated_traces:
        if current in color_dict:  # Use color from the color dictionary
            ax2.plot(tInj_np, eInjVec_np, label=f'{current*1000} pA', linestyle='-')
            ax2.set_xlabel('Time (ms)')
            ax2.set_ylabel('Membrane Potential (mV)')
            ax2.set_title('Simulated Traces')
            ax2.legend(loc='upper right', fontsize='small', frameon=False)
            ax2.set_ylim(-130, 40)
            ax3.plot(tInj_np, eInjVec_np, label=f'{current*1000} pA', linestyle='-')

    # Set labels and title
    ax3.set_xlabel('Time (ms)')
    ax3.set_ylabel('Membrane Potential (mV)')
    ax3.set_title('Overlayed Traces')

    plt.show()

    # Execution time
    end = clock.time()
    print("Execution time of the program is- ", end-start)

    # Return data for further analysis if needed
    return simulated_traces, normalized_hyperpol_current, time_points

def single_comp_model_VC_nachbac(initial_mV, voltages=None, continueRun=None, inject_time=None, injDur=None, delay = None):
    """Function is a simplified voltage clamp protocal to generate an I-V curve for the NaChBac channel 
        comparing back to the I-V curve from Strege et al 2023 (https://doi.org/10.7554/eLife.79271)
  
    Parameters:
        initial_mV: resting membrane potential to initialize the model from
        continueRun: How long the simulation will go for 
        injDur: Duration of the injection time 
        inject_time: preclamp injection time
        delay: Delay of time before current injections (ms)
        

        """
    # Define the soma and its properties
    soma = h.Section(name="soma")
    soma.L = 6.98
    soma.diam = 5.50
    # soma.insert('pas')
    soma.insert('na_bac_strege')
    soma.ena = 60
    # for seg in soma:
        # seg.pas.g = 0.0003796
        # seg.pas.e = -65

    voltage_steps = voltages  # From -100 mV to +30 mV in 10 mV increments

    # Simulation parameters
    delay = delay  # ms, delay before the clamp starts
    initial_mV = initial_mV  # Initial resting potential
    tstop = continueRun  # ms, total simulation time
    simulated_traces = []

    inject_time = inject_time

    # Set up the voltage clamp
    vclamp = h.VClamp(soma(0.5))
    vclamp.dur[0] = inject_time  
    vclamp.amp[0] = -100  # Pre-clamp holding potential


    vclamp.dur[1] = injDur  # Injection duration for the voltage step
    vclamp.amp[2] = -100  
    vclamp.dur[2] = 100  # Post-clamp duration
    h.tstop = tstop

    # Vectors for recording
    tVec = h.Vector().record(h._ref_t)  # Time
    iVec = h.Vector().record(vclamp._ref_i)  # Clamp current
    vVec = h.Vector().record(soma(0.5)._ref_v)  # Voltage at soma

    # Run simulation for each voltage step
    step_voltage = []
    plt.figure()
    for step in voltage_steps:
        vclamp.amp[1] = step  # Set the voltage for the step
        
        # Initialize and run the simulation
        h.finitialize(initial_mV)
        h.continuerun(tstop)

        # Retrieve and store simulation data
        t = np.array(tVec)
        i = np.array(iVec)
        v = np.array(vVec)
        # v_e = np.array(eVec)
        simulated_traces.append({'time': t, 'current': i, 'voltage_soma': v, 'step': step})

        # Retrieve the protocol (durations and amplitudes for this step)
        durations = [vclamp.dur[i] for i in range(len(vclamp.dur))]  # Use len() here
        amplitudes = [vclamp.amp[i] for i in range(len(vclamp.amp))]  # Use len() here

        # Generate time points and voltage values for the protocol
        time_points = [0]  # Start at 0 ms
        voltage_trace = []  # Voltage values over time

        for dur, amp in zip(durations, amplitudes):
            time_points.append(time_points[-1] + dur)  # Add duration to previous time point
            voltage_trace.append(amp)  # Add the corresponding amplitude for this step

        # Create step-like data for plotting
        step_time = []
        step_voltage = []

        for i in range(len(voltage_trace)):
            step_time.extend([time_points[i], time_points[i + 1]])
            step_voltage.extend([voltage_trace[i], voltage_trace[i]])

    # Plot the voltage step protocol for this simulation
        plt.plot(step_time, step_voltage, label=f"Step: {step} mV", drawstyle='steps-post')
        print (f"done {step} mV voltage step")
        
    # Add plot labels and legend
    plt.xlabel("Time (ms)")
    plt.ylabel("Voltage (mV)")
    plt.title("Voltage Clamp Protocol for All Steps")
    plt.grid()
    plt.legend()
    print("done")
    plt.show()
    print("done")

    # Example: Plot results for one voltage step

    # plt.figure(figsize=(10, 5))
    for trace in simulated_traces:
        plt.plot(trace['time'], trace['current']*100, label=f'{trace["step"]} mV')
    plt.xlim(inject_time+injDur-50, continueRun-130)
    plt.ylim(-400, 400)
    plt.xlabel('Time (ms)',fontsize=12)
    plt.ylabel('Clamp Current (pA)',fontsize=12)
    plt.title('Voltage Clamp Experiment',fontsize=12)
    plt.legend(
    loc='lower center',
    bbox_to_anchor=(0.45, -0.41),  # Center horizontally, position below
    ncol=5  # Arrange legend in 3 columns (adjust as needed)
    )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.26)
    plt.show()

    plot_iv_curve(simulated_traces, injDur, delay, inject_time)


def activateCurrentInjectionSOMA_gating_kinetics(initial_mV, sizSection, somaSection, electrodeSec, continueRun=150, current=None, injDur=None, delay= None, type = None, plots = False):
    start = clock.time()
    h.dt = 0.1
    coreneuron.enable = True
    coreneuron.verbose = 0
    coreneuron.model_stats = True
    coreneuron.num_gpus = 1
    h.v_init = initial_mV   # Resting membrane potential 

    if type == "nachbac":
        na_bac_gate_siz = h.Vector()
        na_bac_gate_soma = h.Vector()
        na_bac_h_siz = h.Vector()
        na_bac_m_siz = h.Vector()
        na_bac_m_siz.record(sizSection(0.05)._ref_m_na_bac_strege)
        na_bac_h_siz.record(sizSection(0.05)._ref_h_na_bac_strege)
        na_bac_gate_siz.record(sizSection(0.05)._ref_ina_na_bac_strege)
        na_bac_gate_soma.record(somaSection(0.05)._ref_ina2_na_bac_strege)
        na_bac_current_siz = h.Vector()
        na_bac_current_siz.record(sizSection(0.05)._ref_ina_na_bac_strege)

    elif type == "wt":
        pass
    
    t_vec = h.Vector()
    t_vec.record(h._ref_t)

    ks_gunay_m_siz = h.Vector()
    nap_m_siz = h.Vector()
    nat_m_siz = h.Vector()
    nat_h_siz = h.Vector()

    ks_gunay_m_siz.record(sizSection(0.05)._ref_m_ks_gunay)
    nap_m_siz.record(sizSection(0.05)._ref_m_nap)
    nat_m_siz.record(sizSection(0.05)._ref_m_nat)
    nat_h_siz.record(sizSection(0.05)._ref_h_nat)
        
    # Measuring current of channel and probability of opening
    nap_gate_siz = h.Vector()
    nat_gate_siz = h.Vector()
    ks_gunay_gate_siz = h.Vector()



    nap_gate_siz.record(sizSection(0.05)._ref_ina2nap_nap)
    nat_gate_siz.record(sizSection(0.05)._ref_ina2nat_nat)
    ks_gunay_gate_siz.record(sizSection(0.05)._ref_ik2_ks_gunay)

    #Measuring current of individual channels
    nap_current_siz = h.Vector()
    nachbac_current_siz = h.Vector()
    nat_current_siz = h.Vector()
    ks_gunay_current_siz = h.Vector()
    
    nap_current_siz.record(sizSection(0.05)._ref_ina_nap)
    nat_current_siz.record(sizSection(0.05)._ref_ina_nat)
    ks_gunay_current_siz.record(sizSection(0.05)._ref_ik_ks_gunay)


    stimobj = h.IClamp(somaSection(0.5))
    stimobj.delay = delay
    stimobj.dur = injDur
    stimobj.amp = current
    stimobj.i = 0

    SIZsec = h.Vector()
    SIZsec.record(sizSection(0.02)._ref_v)

    Somasec = h.Vector().indgen()
    Somasec.record(somaSection(0.5)._ref_v)
    Electrodesec = h.Vector()
    Electrodesec.record(electrodeSec(0.5)._ref_v)

    h.finitialize(initial_mV * mV)

    # Run the simulation
    h.tstop = continueRun  # ms
    cnargs = coreneuron.nrncore_arg(h.tstop)
    run()

    # Arguments for CoreNEURON run

    # Re-initialize and run CoreNEURON
    # h.stdinit()
    pc.barrier()
    pc.nrncore_run(cnargs, 1)

    end = clock.time()
    print("Execution time of the program is- ", end-start)


    # Convert NEURON vectors to NumPy arrays
    t_numpy = np.array(t_vec)
    Somasec_numpy = np.array(Somasec)

    # Save data to CSV
    # filename = f'current_{type}_{current:.2f}.csv'
    # with open(filename, 'w', newline='') as csvfile:
    #     csvwriter = csv.writer(csvfile)
    #     csvwriter.writerow(['Time (ms)', 'Somasec (mV)'])
    #     for t_val, soma_val in zip(t_numpy, Somasec_numpy):
    #         csvwriter.writerow([t_val, soma_val])



    if plots == True:
        fig, axs = plt.subplots(1, 3, figsize=(12, 4))  # 1 row, 3 columns of subplots
        # Plot 1
        axs[0].plot(t_vec, SIZsec, 'r')
        axs[0].plot(t_vec, Somasec, 'g')
        axs[0].plot(t_vec, Electrodesec, 'b')
        axs[0].axhline(y=0, color='black', linestyle='-')
        axs[0].axhline(y=initial_mV + 10, color='black', linestyle='--')
        axs[0].set_xlabel('Time (ms)')
        axs[0].set_ylabel('Voltage (mV)')
        axs[0].legend(['SIZ', 'Soma', "ElectrodeSec", '0 mV', '10mV from resting'])
        axs[0].set_title(f'Depolarization of KC in response to SOMA {current} pA current injection')

        # Plot 2
        
        # axs[3].plot(t_vec, ks_gunay_gate_siz, label='ks_siz')
        # axs[3].plot(t_vec, nat_gate_siz, label='nat_siz')
        # axs[3].plot(t_vec, nap_gate_siz, label='nap_siz')
        # axs[3].set_xlabel('Time (ms)')
        # axs[3].set_ylabel('Activation')
        # axs[3].set_title('Conductance')
        # axs[3].legend(loc='best')

        # Plot 3
        # axs[3].plot(t_vec, ks_gunay_current_siz, label='ks_siz')
        # axs[3].plot(t_vec, nat_current_siz, label='nat_siz')
        # axs[3].plot(t_vec, nap_current_siz, label='nap_siz')
        # axs[3].plot(t_vec, na_bac_gate_siz, label='nacbac_siz')
        # axs[3].set_xlabel('Time (ms)')
        # axs[3].set_ylabel('Current (mA/cm²)')  # Adjust unit if needed
        # axs[3].set_title('Current of individual channels')
        # axs[3].legend(loc='best')

        # Plot 4
        axs[1].plot(t_vec, nat_m_siz, label='nat_m_siz')
        axs[1].plot(t_vec, nat_h_siz, label='nat_h_siz')
        axs[1].plot(t_vec, nap_m_siz, label='nap_m_siz')
        axs[1].plot(t_vec, ks_gunay_m_siz, label='ks_gunay_m_siz')
        axs[1].set_xlabel('Time (ms)')
        axs[1].set_ylabel('Activation')
        axs[1].set_title('Activation gates')

        axs[2].plot(t_vec, nat_m_siz, label='nat_m_siz')
        axs[2].plot(t_vec, nat_h_siz, label='nat_h_siz')
        axs[2].plot(t_vec, nap_m_siz, label='nap_m_siz')
        axs[2].plot(t_vec, ks_gunay_m_siz, label='ks_gunay_m_siz')
        axs[2].set_xlabel('Time (ms)')
        axs[2].set_ylabel('Activation')
        axs[2].set_title('Activation gates')
        axs[2].set_xlim(200, 500)


        if type == "nachbac":
            # axs[1].plot(t_vec, na_bac_gate_siz, label='na_bac_siz')
            # axs[2].plot(t_vec, na_bac_current_siz, label='na_bac_siz')
            axs[1].plot(t_vec, na_bac_m_siz, label='na_bac_m_siz')
            axs[1].plot(t_vec, na_bac_h_siz, label='na_bac_h_siz')
            axs[1].legend(loc='best')
            axs[2].plot(t_vec, na_bac_m_siz, label='na_bac_m_siz')
            axs[2].plot(t_vec, na_bac_h_siz, label='na_bac_h_siz')
            axs[2].legend(loc='best')
        elif type == "wt":
            axs[1].legend(loc='best')
            # axs[2].legend(loc='best')
        
        
        # Adjust layout
        plt.tight_layout()
        plt.show()
    #For individual plots

    # plt.figure(figsize=(12, 8))
    # plt.plot(t_vec, SIZsec, 'r')
    # plt.plot(t_vec, Somasec, 'g')
    # plt.axhline(y=0, color='black', linestyle='-')
    # plt.axhline(y=initial_mV+10, color='black', linestyle='--')
    # plt.xlabel('Time (ms)')
    # plt.ylabel('Voltage (mV)')
    # plt.plot(t_vec, Electrodesec, 'b')
    # plt.legend(['SIZ', 'Soma', '0 mV', '10mV from resting', "ElectrodeSec"])
    # maxSIZ_depol = SIZsec.max()
    # plt.title('Depolarization of KC in response to SOMA {} pA current injection, SIZ depolarized to {} mV, SIZ depolarized to {} mV above resting'.format(current, round(maxSIZ_depol, 3), round(SIZsec.max() - initial_mV, 3)))

    # plt.figure(figsize=(12, 8))
    # plt.plot(t_vec, Electrodesec, 'red')
    # plt.plot(t_vec, SIZsec, 'green')
    # plt.plot(t_vec, SIZsec2, 'blue')
    # plt.plot(t_vec, SIZsec3, 'black')
    # plt.xlabel('Time (ms)')
    # plt.ylabel('Voltage (mV)')
    # plt.legend(['Electrode', 'SIZ(0.12)','SIZ(0.5)', 'SIZ(0.99)'])
    # plt.title('Recordings from different SIZ locations + Electrode')

    # plt.figure(figsize=(12, 8))
    # plt.plot(t_vec, Electrodesec, 'black')
    # plt.plot(t_vec, axon1, 'red')
    # plt.plot(t_vec, axon2, 'orange')
    # plt.plot(t_vec, axon3, 'yellow')
    # plt.plot(t_vec, axon4, 'green')
    # plt.plot(t_vec, axon5, 'blue')
    # plt.plot(t_vec, axon6, 'violet')
    # plt.plot(t_vec, axon7, 'cyan')
    # plt.xlabel('Time (ms)')
    # plt.ylabel('Voltage (mV)')
    # plt.legend(['Electrode', 'Axon5','Axon15', 'Axon25', 'Axon35', 'Axon45', 'Axon55', 'Axon125'])
    # plt.title('Recordings from different SIZ locations + Electrode')
 
    # plt.figure(figsize=(12, 8))
    # plt.plot(t_vec, na_bac_gate_siz, label='na_bac_siz')
    # plt.plot(t_vec, na_bac_gate_soma, label='na_bac_soma')
    # plt.plot(t_vec, ks_gunay_gate_siz, label='ks_siz')
    # plt.plot(t_vec, nat_gate_siz, label='nat_siz')
    # plt.plot(t_vec, nap_gate_siz, label='nap_siz')
    # plt.xlabel('Time (ms)')
    # plt.ylabel('Activation')
    # plt.title('Activation gates')
    # plt.legend(loc='best')  # Automatically finds the best position for the legend
   
    # plt.figure(figsize=(12, 8))
    # plt.plot(t_vec, na_bac_current_siz, label='na_bac_siz')
    # plt.plot(t_vec, na_bac_current_soma, label='na_bac_soma')
    # plt.plot(t_vec, ks_gunay_current_siz, label='ks_siz')
    # plt.plot(t_vec, nat_current_siz, label='nat_siz')
    # plt.plot(t_vec, nap_current_siz, label='nap_siz')
    # plt.xlabel('Time (ms)')
    # plt.ylabel('Current (mA/cm2) <- need to check official units')
    # plt.title('Current of individual channels')
    # plt.legend(loc='best')  # Automatically finds the best position for the legend

    # plt.figure(figsize=(12, 8))
    # plt.plot(t_vec, na_bac_m_soma, label='na_bac_m_soma')
    # plt.plot(t_vec, na_bac_h_soma, label='na_bac_h_soma')
    # plt.plot(t_vec, na_bac_m_siz, label='na_bac_m_siz')
    # plt.plot(t_vec, na_bac_h_siz, label='na_bac_h_siz')
    # plt.plot(t_vec, nat_m_siz, label='nat_m_siz')
    # plt.plot(t_vec, nat_h_siz, label='nat_h_siz')
    # plt.plot(t_vec, nap_m_siz, label='nap_m_siz')
    # plt.plot(t_vec, ks_gunay_m_siz, label='ks_gunay_m_siz')
    # plt.xlabel('Time (ms)')
    # plt.ylabel('Activation')
    # plt.title('Activation gates')
    # plt.legend(loc='best')  # Automatically finds the best position for the legend
   
    # plt.figure(figsize=(12, 8))
    # plt.plot(t_vec, na_bac_m_soma, label='na_bac_m_soma')
    # plt.plot(t_vec, na_bac_h_soma, label='na_bac_h_soma')
    # plt.plot(t_vec, na_bac_m_siz, label='na_bac_m_siz')
    # plt.plot(t_vec, na_bac_h_siz, label='na_bac_h_siz')
    # plt.xlabel('Time (ms)')
    # plt.ylabel('Activation')
    # plt.title('Activation gates')
    # plt.legend(loc='best')


#####################

#Stastitical testing

def run_stats(column=None):
    """
    Functions to run stats, checking normality, via shapiro-wilk, equality of variance via levenes test, and t-test based on levenes results
    Parameter:
        column: which column within the parameters list to compare across model types ex. nat_siz, nap_siz, nat_axon etc
    """
    if column is None:
        raise ValueError("Please provide a column name to analyze.")

    current_dir = os.getcwd()
    file_name = 'Simulation Passive and Active properties.xlsx'  # Adjust as needed
    file_path = os.path.join(current_dir, file_name)

    # Read the Excel file
    data = pd.read_excel(file_path)

    # Strip column names of any extra spaces
    data.columns = data.columns.str.strip()
    
    # Filter data for nachbac and WT only
    subset_data = data[data['Type_simulation'].isin(['nachbac', 'WT'])]

    # Check normality using Shapiro-Wilk test
    stat_nachbac, p_nachbac = stats.shapiro(subset_data[subset_data['Type_simulation'] == 'nachbac'][column])
    stat_wt, p_wt = stats.shapiro(subset_data[subset_data['Type_simulation'] == 'WT'][column])

    print(f"Shapiro-Wilk test for nachbac: p-value = {p_nachbac:.4f}")
    print(f"Shapiro-Wilk test for WT: p-value = {p_wt:.4f}")

    # Determine if data is normally distributed
    normal_nachbac = p_nachbac > 0.05
    normal_wt = p_wt > 0.05

    # Perform statistical test based on normality
    if normal_nachbac and normal_wt:
        # Check for equal variance using Levene's test
        levene_stat, levene_p = stats.levene(
            subset_data[subset_data['Type_simulation'] == 'nachbac'][column],
            subset_data[subset_data['Type_simulation'] == 'WT'][column]
        )
        print(f"Levene's test for equal variances: p-value = {levene_p:.4f}")

        if levene_p > 0.05:
            # Perform independent t-test (equal variances assumed)
            t_stat, p_ttest = stats.ttest_ind(
                subset_data[subset_data['Type_simulation'] == 'nachbac'][column],
                subset_data[subset_data['Type_simulation'] == 'WT'][column],
                equal_var=True
            )
            print(f"Independent t-test: p-value = {p_ttest:.4f}")
        else:
            # Perform Welch’s t-test (unequal variances)
            t_stat, p_ttest = stats.ttest_ind(
                subset_data[subset_data['Type_simulation'] == 'nachbac'][column],
                subset_data[subset_data['Type_simulation'] == 'WT'][column],
                equal_var=False
            )
            print(f"Welch's t-test (unequal variances): p-value = {p_ttest:.4f}")
    else:
        # Perform Mann-Whitney U test if normality is violated
        u_stat, p_mannwhitney = stats.mannwhitneyu(
            subset_data[subset_data['Type_simulation'] == 'nachbac'][column],
            subset_data[subset_data['Type_simulation'] == 'WT'][column],
            alternative='two-sided'
        )
        print(f"Mann-Whitney U test: p-value = {p_mannwhitney:.4f}")

    # Boxplot with mean values overlaid
    plt.figure(figsize=(8, 5))
    
    # Boxplot
    sns.boxplot(data=subset_data, x='Type_simulation', y=column, boxprops=dict(facecolor='none'))

    # Stripplot for individual data points
    sns.stripplot(data=subset_data, x='Type_simulation', y=column, color='blue', alpha=0.6, jitter=True)

    # Overlay mean values in red
    sns.pointplot(
        data=subset_data, x='Type_simulation', y=column, 
        estimator=np.mean, color='red', join=False, markers='o', scale=1.2
    )

    plt.title(f'Comparison of {column} between nachbac and WT')
    plt.xlabel('Type Simulation')
    plt.ylabel(column)
    plt.tight_layout()
    plt.show()

#####################

def main():
    morph_file = ("720575940606954507_skeletonized_um_model.swc") #morphology swc model

    #WT passive property sets:
    #wt fly 8
    # raVal = 110
    # gleakVal = 0.000110
    # cmVal = 2
    # erev = -78.25

    #wt fly 9
    # raVal = 100
    # gleakVal = 0.0001115
    # cmVal = 2
    # erev = -71.25
    
    #wt fly 10
    # raVal = 128
    # gleakVal = 0.000110
    # cmVal = 2
    # erev = -80.25

    #NaChBac fly 3 
    # raVal = 75
    # gleakVal = 0.000038755
    # cmVal = 1
    # erev = -88

    #NaChBac fly 6
    # raVal = 70
    # gleakVal = 0.0023215
    # cmVal = 2
    # erev = -62.5

    #NaChBac fly 14 
    raVal = 125
    gleakVal = 0.0000705
    cmVal = 2
    erev = -85.75 

    cell, allSections_py, allSections_nrn, somaSection, sizSection, erev, axonList, tetherList, dendList, electrodeSec = initializeModel(morph_file, resting_ev= -60.25, Channels=True,  capacitance = cmVal, axial_resistivity = raVal, gleak_conductance= gleakVal, erev = erev)


    elec_raVal = 205.513  
    elec_cmVal = 6.4

    sealCon_8GOhm = 0.0003978 
    elec_gleakVal = sealCon_8GOhm

    change_Ra(ra=raVal, electrodeSec=electrodeSec, electrodeVal = elec_raVal)
    change_gLeak(gleak=gleakVal, erev=erev, electrodeSec=electrodeSec, electrodeVal = elec_gleakVal)
    change_memCap(memcap=cmVal, electrodeSec=electrodeSec, electrodeVal = elec_cmVal)


    ### Running the simulations with injecting current at the SOMA
    #To change the active properties do that in the initializeModel function.

    #Fitting passive properties of multiple current injections
    resting_ev = -68 #The initializing membrane potential
    # currents_list = [-0.04, -0.03, -0.02, -0.01, 0]
    # simulated_traces, time_points_np, current_export = passive_fitting(resting_ev, electrodeSec, somaSection, currents=currents_list, continueRun=2000, injDur=1000, delay = 231.4, type = "wt", flynum = 9)
    
    # #RMSE for beginning and end of current injections to fit time constants, calculates RMSE for the beginning and ending of the current injection matching the time constants. 
    # Calculate_RMSE(start_index=11000, end_index=16000, time_data=time_points_np, current_data=current_export, simulated_traces=simulated_traces, trace_index=3)
    # Calculate_RMSE(start_index=57000, end_index=70000, time_data=time_points_np, current_data=current_export, simulated_traces=simulated_traces, trace_index=3)
    
    #####################

    #Simulations and plotting

    #Simulation for Fig 7, panel C-G

    currents_list = [-0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, .04, .05, .06]
    runSimAndNormalize(initial_mV= resting_ev , electrodeSec=electrodeSec, somaSection=somaSection, currents=currents_list, continueRun=2000, injDur=1000, delay=231.4, type='nachbac', flynum = 14, current_start = -40, current_end = 60)
    
    #Statistics and plotting for Fig 7, panel H
    # run_stats(column = 'nat_siz')

    #Simulations for figure 7, panel I
    # currents_list = [-0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10] 
    # Soma_current_injections(resting_ev, electrodeSec, somaSection, currents=currents_list, continueRun=2000, injDur=1000, delay = 231.4, type = "nachbac",  manipulation = True, flynum= 6)
    
    #Plotting of figure 7, panel I
    # plot_spike_counts_vs_currents("wt_simulation_results.csv",  "simulation_results_nachbac_final.csv")
    
    #Simluation for figure 7, panel J
    # currents_list = [-0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, .04, .05, .06, 0.07, 0.08, 0.09, .1]
    # Varying_nachbac_soma_current_injections(resting_ev, electrodeSec, somaSection, currents=currents_list, continueRun=2000, injDur=1000, delay = 231.4, type = "nachbac")

    #Plotting for figure 7 panel J
    # plot_spike_count_vs_current('WT_simulation_with_varying_nachbac_cond.csv')

    #Plotting of figure # panel ...

    # activateCurrentInjectionSOMA_gating_kinetics(resting_ev, sizSection, somaSection, electrodeSec, continueRun=2000, current=0.02, injDur=1000, delay= 231.4, type = "nachbac", plots = True)
    

#Runs the script, uncomment simulations as needed 
main()