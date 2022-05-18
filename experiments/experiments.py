import sys

sys.path.append("../")

from cnot import standard_cnot
from qiskit import Aer, execute
from qiskit.circuit import QuantumCircuit, QuantumRegister
from qiskit.providers.aer import AerSimulator
from qiskit.providers.aer.noise import NoiseModel
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

    
def counts_experiment(circuit: QuantumCircuit, shots: int, noise_model: NoiseModel = None) -> dict:
    """Executes a sampling experiment.

    Args:
        circuit (QuantumCircuit): A quantum circuit.
        shots (int): The number of shots of the experiment.
        noise_model (NoiseModel, optional): The noise model of the simulation. Defaults to None.

    Returns:
        dict: The dictionary of counts.
    """
    if isinstance(circuit, QuantumCircuit) and isinstance(shots, int) and (isinstance(noise_model, NoiseModel) or noise_model is None):
        
        counts = execute(experiments = circuit, backend = AER_SIM, shots = shots, noise_model = noise_model).result().get_counts()
        return counts
        
    else:
        
        raise TypeError("The inputs are not a QuantumCircuit, an int, and a NoiseModel.")


def cnot_truth_table_experiment(circuit: QuantumCircuit) -> dict:
    """Tests if a given circuit represents a CNOT gate.

    Args:
        circuit (QuantumCircuit): A circuit of a candidate a CNOT gate.

    Raises:
        TypeError: If the input is not a QuantumCircuit.

    Returns:
        dict: A dictionary with the results of the experiment.
    """
    if isinstance(circuit, QuantumCircuit):
        
        inputs = ['00', '01', '10', '11']
        experiments_results_dict = {}
        
        for input in inputs:
            
            if input == '00':
                
                counts = counts_experiment(circuit = circuit, shots = 1)
                result = '00' in counts.keys()
                experiments_results_dict[input+" -> 00"] = result
            elif input == '01':
                
                qc = QuantumCircuit(2,2)
                qc.x(qubit = 0)
                qc = qc.compose(circuit)
                counts = counts_experiment(circuit = qc, shots = 1)
                result = '11' in counts.keys()
                experiments_results_dict[input+" -> 11"] = result
            elif input == '10':
                
                qc = QuantumCircuit(2,2)
                qc.x(qubit = 1)
                qc = qc.compose(circuit)
                counts = counts_experiment(circuit = qc, shots = 1)
                result = '10' in counts.keys()
                experiments_results_dict[input+" -> 10"] = result
            else:
                
                qc = QuantumCircuit(2,2)
                qc.x(qubit = [0,1])
                qc = qc.compose(circuit)
                counts = counts_experiment(circuit = qc, shots = 1)
                result = '01' in counts.keys()
                experiments_results_dict[input+" -> 01"] = result
                
        return experiments_results_dict      
                
    else:
        
        raise TypeError("The input is not a QuantumCircuit!")

    
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
                    
                    if keys.split(" ")[-1] == logical_bit*code_qubits_register.size:
                        
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