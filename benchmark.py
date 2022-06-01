import datetime as dt
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings("ignore")

from cnot import controlled_rx_cnot, floating_gate_cnot, hadamard_cz_cnot, molmer_sorensen_cnot, standard_cnot, sqrt_xx_cnot, sqrt_swap_cnot
from experiments import counts_experiment, process_result, process_infidelity_experiment, avg_gate_infidelity_experiment, cnot_truth_table_experiment
from qiskit.providers.aer.noise import depolarizing_error, pauli_error
from plot import plot_experiments_results
from qec import create_3_repetition_code_circuit

date_and_time = str(dt.datetime.now()).replace(":","-").split(".")[0].replace(" ", "_")
cnots_list = ["standard", "hadamard_cz", "sqrt_xx", "floating_gate", "controlled_rx", "molmer_sorensen", "sqrt_swap"]

columns_truth_table_experiment = ["type", "00 -> 00", "01 -> 11", "10 -> 10", "11 -> 01"]
data_truth_table_experiment = []

for cnot in cnots_list:
    
    aux_list = []
    aux_list.append(cnot)
    if cnot == "standard":
        
        circuit = standard_cnot(measurements = True)
        result = cnot_truth_table_experiment(circuit = circuit)
        aux_list2 = [value for _, value in result.items()]
        aux_list.extend(aux_list2)
    elif cnot == "hadamard_cz":
        
        circuit = hadamard_cz_cnot(measurements = True)
        result = cnot_truth_table_experiment(circuit = circuit)
        aux_list2 = [value for _, value in result.items()]
        aux_list.extend(aux_list2)
    elif cnot == "sqrt_xx":
        
        circuit = sqrt_xx_cnot(measurements = True)
        result = cnot_truth_table_experiment(circuit = circuit)
        aux_list2 = [value for _, value in result.items()]
        aux_list.extend(aux_list2)
    elif cnot == "floating_gate":
        
        circuit = floating_gate_cnot(measurements = True)
        result = cnot_truth_table_experiment(circuit = circuit)
        aux_list2 = [value for _, value in result.items()]
        aux_list.extend(aux_list2)
    elif cnot == "controlled_rx":
        
        circuit = controlled_rx_cnot(measurements = True)
        result = cnot_truth_table_experiment(circuit = circuit)
        aux_list2 = [value for _, value in result.items()]
        aux_list.extend(aux_list2)
    elif cnot == "molmer_sorensen":
        
        circuit = molmer_sorensen_cnot(measurements = True)
        result = cnot_truth_table_experiment(circuit = circuit)
        aux_list2 = [value for _, value in result.items()]
        aux_list.extend(aux_list2)
    else:
        
        circuit = sqrt_swap_cnot(measurements = True)
        result = cnot_truth_table_experiment(circuit = circuit)
        aux_list2 = [value for _, value in result.items()]
        aux_list.extend(aux_list2)
        
    data_truth_table_experiment.append(aux_list)
    
df_data_truth_table_experiment = pd.DataFrame(data = data_truth_table_experiment, columns = columns_truth_table_experiment)
df_data_truth_table_experiment.to_csv("results/checking_truth_table_experiment.csv")

columns = ["Error probability", "00 Counts", "01 Counts", "10 Counts", "11 Counts", "p_00", "p_not_00", "Process infidelity", "Average gate infidelity"]
columns_qec = ["Error probability", "Process infidelity", "Average gate infidelity", "p_success", "p_logical_error"]
type_of_errors_list = ["bitflip", "depolarizing", "bitflip_and_coherent_errors", "depolarizing_and_coherent_errors"]

errors_list = sorted([0.0, 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001, 0.00075, 0.0005, 0.00025, 0.0001])

controlled_rx_params = [np.pi + ((-1)**np.random.randint(0,2))*np.random.random()]
molmer_sorensen_params = [np.pi/2 + ((-1)**np.random.randint(0,2))*np.random.random(), np.pi/2 + ((-1)**np.random.randint(0,2))*np.random.random(), -np.pi/2 + ((-1)**np.random.randint(0,2))*np.random.random(), -np.pi/2 + ((-1)**np.random.randint(0,2))*np.random.random(), -np.pi/2 + ((-1)**np.random.randint(0,2))*np.random.random()]
sqrt_swap_params = [np.pi/2 + ((-1)**np.random.randint(0,2))*np.random.random(), -np.pi/2 + ((-1)**np.random.randint(0,2))*np.random.random(), -np.pi/2 + ((-1)**np.random.randint(0,2))*np.random.random(), -np.pi/2 + ((-1)**np.random.randint(0,2))*np.random.random()]

params_list = [controlled_rx_params, molmer_sorensen_params, sqrt_swap_params]

df_params = pd.DataFrame(data = params_list)
df_params.to_csv(f"results/parameters_used_in_coherent_errors_{date_and_time}.csv")

for cnot in cnots_list:
    
    for error in type_of_errors_list:
        
        data_cnot_infidelity = []
        data_qec_cnot_infidelity = []
    
        for p_error in errors_list:
            
            if error == "bitflip" or error == "bitflip_and_coherent_errors":
                
                bitflip = pauli_error([('X', p_error), ('I', 1 - p_error)])
                cx_error = bitflip.tensor(bitflip)          
            else:
                
                cx_error = depolarizing_error(param = p_error, num_qubits = 2)

            if cnot == "standard":
                
                if error == "bitflip_and_coherent_errors" or error == "depolarizing_and_coherent_errors":
                    
                    continue
                else:
                    
                    qc_with_meas = standard_cnot(measurements = True, noisy = cx_error)
                    qc_without_meas = standard_cnot(noisy = cx_error)
            elif cnot == "hadamard_cz":
                
                if error == "bitflip_and_coherent_errors" or error == "depolarizing_and_coherent_errors":
                    
                    continue
                else:
                    
                    qc_with_meas = hadamard_cz_cnot(measurements = True, noisy = cx_error)
                    qc_without_meas = hadamard_cz_cnot(noisy = cx_error)
            elif cnot == "sqrt_xx":
                
                if error == "bitflip_and_coherent_errors" or error == "depolarizing_and_coherent_errors":
                    
                    continue
                else:
                    
                    qc_with_meas = sqrt_xx_cnot(measurements = True, noisy = cx_error)
                    qc_without_meas = sqrt_xx_cnot(noisy = cx_error)
            elif cnot == "floating_gate":
                
                if error == "bitflip_and_coherent_errors" or error == "depolarizing_and_coherent_errors":
                    
                    continue
                else:
                    
                    qc_with_meas = floating_gate_cnot(measurements = True, noisy = cx_error)
                    qc_without_meas = floating_gate_cnot(noisy = cx_error) 
                
            elif cnot == "controlled_rx":
                
                if error == "bitflip" or error == "depolarizing":
                    
                    qc_with_meas = controlled_rx_cnot(measurements = True, noisy = cx_error)
                    qc_without_meas = controlled_rx_cnot(noisy = cx_error)
                else:
                    
                    qc_with_meas = controlled_rx_cnot(measurements = True, noisy = cx_error, params = params_list[0])
                    qc_without_meas = controlled_rx_cnot(noisy = cx_error, params = params_list[0])
            elif cnot == "molmer_sorensen":
                
                if error == "bitflip" or error == "depolarizing":
                    
                    qc_with_meas = molmer_sorensen_cnot(measurements = True, noisy = cx_error)
                    qc_without_meas = molmer_sorensen_cnot(noisy = cx_error)
                else:
                    
                    qc_with_meas = molmer_sorensen_cnot(measurements = True, noisy = cx_error, params = params_list[1])
                    qc_without_meas = molmer_sorensen_cnot(noisy = cx_error, params = params_list[1])
            else:
                
                if error == "bitflip" or error == "depolarizing":
                    
                    qc_with_meas = sqrt_swap_cnot(measurements = True, noisy = cx_error)
                    qc_without_meas = sqrt_swap_cnot(noisy = cx_error)
                else:
                    
                    qc_with_meas = sqrt_swap_cnot(measurements = True, noisy = cx_error, params = params_list[2])
                    qc_without_meas = sqrt_swap_cnot(noisy = cx_error, params = params_list[2])
                    
            qc_code = create_3_repetition_code_circuit(env_param = np.pi/2, noisy_cnot = qc_without_meas)
            
            counts = counts_experiment(circuit = qc_with_meas, shots = 8192)
            counts_code = counts_experiment(circuit = qc_code, shots = 8192)
            processed_counts = process_result(code_qubits_register = qc_with_meas.qregs[0], logical_bit = "0", counts = counts, num_shots = 8192)
            processed_counts_code = process_result(code_qubits_register = qc_code.qregs[1], logical_bit = "0", counts = counts_code, num_shots = 8192, has_ancilla = True)
            process_infidelity = process_infidelity_experiment(circuit = qc_without_meas)
            avg_gate_infidelity = avg_gate_infidelity_experiment(circuit = qc_without_meas)
                    
            if '00' in counts.keys():
                        
                counts_00 = counts['00']
            else:
                        
                counts_00 = 0
            if '01' in counts.keys():
                        
                counts_01 = counts['01']
            else:
                        
                counts_01 = 0
            if '10' in counts.keys():
                        
                counts_10 = counts['10']
            else:
                        
                counts_10 = 0
            if '11' in counts.keys():
                        
                counts_11 = counts['11']
            else:
                        
                counts_11 = 0
                        
            data_cnot_infidelity.append([p_error, counts_00, counts_01, counts_10, counts_11, processed_counts[0], processed_counts[1], process_infidelity, avg_gate_infidelity])
            data_qec_cnot_infidelity.append([p_error, process_infidelity, avg_gate_infidelity, processed_counts_code[0], processed_counts_code[1]])
        
        if len(data_cnot_infidelity) != 0 and len(data_qec_cnot_infidelity) != 0:
                        
            df_data_cnot_infidelity = pd.DataFrame(data = data_cnot_infidelity, columns = columns)
            df_data_cnot_infidelity.to_csv(f"results/{error}_{cnot}_cnot_{date_and_time}.csv")
            df_data_qec_cnot_infidelity = pd.DataFrame(data = data_qec_cnot_infidelity, columns = columns_qec)
            df_data_qec_cnot_infidelity.to_csv(f"results/qec_313_{error}_{cnot}_cnot_{date_and_time}.csv")
            
            plot_experiments_results(data = df_data_cnot_infidelity, type_data = "infidelity", save_fig = True, fig_name = f"{error}_{cnot}_cnot_{date_and_time}")
            plot_experiments_results(data = df_data_qec_cnot_infidelity, type_data = "qec", save_fig = True, fig_name = f"qec_313_{error}_{cnot}_cnot_{date_and_time}")
        