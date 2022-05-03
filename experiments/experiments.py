import sys

sys.path.append("../")

import pandas as pd

from cnot import standard_cnot
from qiskit import Aer, execute
from qiskit.circuit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.providers.aer import AerSimulator
from qiskit.providers.aer.noise import depolarizing_error, pauli_error
from qiskit.quantum_info import average_gate_fidelity, Kraus, process_fidelity

from typing import Optional

QASM_SIM = Aer.get_backend('qasm_simulator')
UNITARY_SIM = Aer.get_backend('unitary_simulator')
AER_SIM = AerSimulator(zero_threshold = 1e-5)
STANDARD_CNOT_UNITARY = execute(experiments = standard_cnot(), backend = UNITARY_SIM).result().get_unitary()

def avg_gate_infidelity_experiment(circuit: QuantumCircuit) -> float:
    """Calculates the average gate infidelity between the quantum circuit
    that represents a CNOT gate and a perfect standard CNOT gate.

    Args:
        circuit (QuantumCircuit): A quantum circuit that represents a CNOT gate.

    Raises:
        TypeError: If the input is not a circuit.

    Returns:
        float: The average gate infidelity between the input circuit and the
        perfect standard CNOT gate.
    """
    if isinstance(circuit, QuantumCircuit):
        
        channel = Kraus(circuit)
        avg_gate_infidelity = 1 - average_gate_fidelity(channel, STANDARD_CNOT_UNITARY)
        
        return avg_gate_infidelity
    else:
        
        raise TypeError("The input is not a QuantumCircuit!")
    
def counts_experiment(circuit: QuantumCircuit, shots: int) -> dict:
    pass
    
def process_infidelity_experiment(circuit: QuantumCircuit) -> float:
    """Calculates the process infidelity between the input quantum circuit and
    the perfect standard CNOT gate.

    Args:
        circuit (QuantumCircuit): A quantum circuit that represents a CNOT gate.

    Raises:
        TypeError: If the input is not a circuit.

    Returns:
        float: The process infidelity between the input and the perfect standard CNOT.
    """
    if isinstance(circuit, QuantumCircuit):
        
        channel = Kraus(circuit)
        infidelity = 1 - process_fidelity(channel, STANDARD_CNOT_UNITARY)
        
        return infidelity
    else:
        
        raise TypeError("The input is not a QuantumCircuit!")
    
def process_result(code_qubits_register: QuantumRegister, logical_bit: str, counts: dict, num_shots: int, has_ancilla: Optional[bool] = False) -> tuple:
    """Processes an experiment result.

    Args:
        code_qubits_register (QuantumRegister): The quantum register used in the quantum circuit.
        logical_bit (str): It tells if we want to analyse 0 or 1 as the right answer.
        counts (dict): A dictionary with the experiment result.
        num_shots (int): The number of shots using in the experiment.
        has_ancilla (Optional[bool], optional): It tells to the function if the counts has ancilla qubits. Defaults to False.. Defaults to False.

    Raises:
        ValueError: If the logical_bit string is not equal to '0' or '1'.
        TypeError: If an argument is passed with an invalid type. 

    Returns:
        tuple: The probabilities of getting the right answer and an error.
    """
    if isinstance(code_qubits_register, QuantumRegister) and isinstance(logical_bit, str) and isinstance(counts, dict) and isinstance(num_shots, int) and isinstance(has_ancilla, bool):
        
        count_sucess = 0
        if logical_bit == '0' or logical_bit == '1':
            
            for keys, value in counts.items():
                
                if has_ancilla:
                    
                    if keys.split(" ")[code_qubits_register.size] == logical_bit*code_qubits_register.size:
                        
                        count_sucess += value
                else:
                    
                    if keys == logical_bit*code_qubits_register.size:
                        
                        count_sucess += value
                        
            prob_sucess = count_sucess/num_shots
            prob_logical_error = 1 - prob_sucess
            probs_tuple = tuple([prob_sucess, prob_logical_error])
            
            return probs_tuple
        else:
            
            raise ValueError("Wrong value of logical! The logical argument must be equal to '0' or '1'.")
    else:
        
        raise TypeError("Invalids arguments! This function only accepts QuantumRegister, QuantumCircuit, str, dict, int and bool types.")

# def test():
    
#     columns = ["Error probability", "00 Counts", "01 Counts", "10 Counts", "11 Counts", "p_00", "p_not_00", "Process infidelity", "Average gate infidelity"]
#     data = []
    
#     for key, value in experiments_counts_dict.items():
        
#         if '00' in value.keys():
            
#             counts_00 = value['00']
#         else:
            
#             counts_00 = 0
#         if '01' in value.keys():
            
#             counts_01 = value['01']
#         else:
            
#             counts_01 = 0
#         if '10' in value.keys():
            
#             counts_10 = value['10']
#         else:
            
#             counts_10 = 0
#         if '11' in value.keys():
            
#             counts_11 = value['11']
#         else:
            
#             counts_11 = 0
            
#         data.append([float(key), counts_00, counts_01, counts_10, counts_11, experiments_])
        
#     data_pd = pd.DataFrame(data = data, columns = columns)