from qiskit.circuit import ClassicalRegister, QuantumCircuit, QuantumRegister

def create_3_repetition_code_circuit(env_param: float, qubit_error_target: int = 1, noisy_cnot: QuantumCircuit = None) -> QuantumCircuit:
    """_summary_

    Args:
        env_param (float): _description_
        qubit_error_target (int, optional): _description_. Defaults to 1.
        noisy_cnot (QuantumCircuit, optional): _description_. Defaults to None.

    Returns:
        QuantumCircuit: _description_
    """
    env_qubit = QuantumRegister(size = 1, name = "environment")
    code_qubits = QuantumRegister(size = 3, name = "code")
    ancillae_qubits = QuantumRegister(size = 3, name = "ancillae")
    bits = ClassicalRegister(size = 3, name = "bits")
    anc_bit0 = ClassicalRegister(size=1, name='anc_bit0')
    anc_bit1 = ClassicalRegister(size=1, name='anc_bit1')
    anc_bit2 = ClassicalRegister(size=1, name='anc_bit2')

    qc = QuantumCircuit(env_qubit, code_qubits, ancillae_qubits, bits, anc_bit0, anc_bit1, anc_bit2)
    
    qc.ry(theta = env_param, qubit = env_qubit)
    
    qc.barrier()
    qc.h(qubit = code_qubits[0])
    qc.cx(control_qubit = code_qubits[0], target_qubit = code_qubits[1])
    qc.cx(control_qubit = code_qubits[0], target_qubit = code_qubits[2])
    qc.barrier()
    
    qc.cx(control_qubit = env_qubit, target_qubit = qubit_error_target)
    
    qc.barrier()
    if noisy_cnot is None:
        
        qc.cx(control_qubit = code_qubits[0], target_qubit = ancillae_qubits[0])
        qc.cx(control_qubit = code_qubits[1], target_qubit = ancillae_qubits[0])
    else:
        
        qc = qc.compose(noisy_cnot, qubits = [1,4])
        qc.cx(control_qubit = code_qubits[1], target_qubit = ancillae_qubits[0])
    qc.barrier()
    
    qc.cx(control_qubit = code_qubits[1], target_qubit = ancillae_qubits[1])
    qc.cx(control_qubit = code_qubits[2], target_qubit = ancillae_qubits[1])
    
    qc.barrier()
    qc.mct(control_qubits = ancillae_qubits[0:2], target_qubit = ancillae_qubits[2])
    qc.cx(control_qubit = ancillae_qubits[2], target_qubit = ancillae_qubits[0])
    qc.cx(control_qubit = ancillae_qubits[2], target_qubit = ancillae_qubits[1])
    qc.barrier()
    
    qc.measure(qubit = ancillae_qubits[0], cbit = anc_bit0)
    qc.measure(qubit = ancillae_qubits[1], cbit = anc_bit1)
    qc.measure(qubit = ancillae_qubits[2], cbit = anc_bit2)
    
    qc.barrier()
    qc.x(qubit = code_qubits[0]).c_if(classical = anc_bit0, val = 1)
    qc.x(qubit = code_qubits[2]).c_if(classical = anc_bit1, val = 1)
    qc.x(qubit = code_qubits[1]).c_if(classical = anc_bit2, val = 1)
    qc.barrier()
    
    qc.cx(control_qubit = code_qubits[0], target_qubit = code_qubits[2])
    qc.cx(control_qubit = code_qubits[0], target_qubit = code_qubits[1])
    qc.h(qubit = code_qubits[0])
    
    qc.barrier()
    qc.measure(qubit = code_qubits, cbit = bits)
    
    return qc


def create_shor_code_circuit() -> QuantumCircuit:
    pass
    