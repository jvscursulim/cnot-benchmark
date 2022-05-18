from qiskit.circuit import ClassicalRegister, QuantumCircuit, QuantumRegister

def create_3_repetition_code_circuit(env_param: float, qubit_error_target: int = 1, noisy_cnot: QuantumCircuit = None) -> QuantumCircuit:
    """Creates the circuit that implements the [[3,1,3]] repetition code.

    Args:
        env_param (float): The angle of the Ry that is applied to the environment qubit.
        qubit_error_target (int, optional): The number of the qubit that will be flipped. Defaults to 1.
        noisy_cnot (QuantumCircuit, optional): A quantum circuit that works like a noisy CNOT gate. Defaults to None.

    Raises:
        TypeError: If the inputs are not equal to: float, int and QuantumCircuit.

    Returns:
        QuantumCircuit: The circuit of [[3,1,3]] repetition code.
    """
    if isinstance(env_param, float) and isinstance(qubit_error_target, int) and (isinstance(noisy_cnot, QuantumCircuit) or noisy_cnot is None):
    
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
    else:
        
        raise TypeError("This function expects as input a float, an int and a QuantumCircuit!")


def create_shor_code_circuit(env_param: float, qubit_error_target: int = 1, noisy_cnot: QuantumCircuit = None) -> QuantumCircuit:
   """Creates a quantum circuit that implements the Shor code [[9,1,3]].

    Args:
        env_param (float): The angle of the Ry that is applied to the environment qubit.
        qubit_error_target (int, optional): The number of the qubit that will be flipped. Defaults to 1.
        noisy_cnot (QuantumCircuit, optional): A quantum circuit that works like a noisy CNOT gate. Defaults to None.

    Returns:
        QuantumCircuit: The quantum circuit of Shor code.
   """  
   env_qubit = QuantumRegister(size = 1, name = "environment")
   code_qubits = QuantumRegister(size = 9, name = "code")
   ancillae_qubits = QuantumRegister(size = 8, name = "ancillae")
   aux_qubits = QuantumRegister(size = 4, name = "auxiliar")
   bits = ClassicalRegister(size = 9, name = "bits")
   classical1 = ClassicalRegister(size = 1, name = "classical1")
   classical2 = ClassicalRegister(size = 1, name = "classical2")
   classical3 = ClassicalRegister(size = 1, name = "classical3")
   classical4 = ClassicalRegister(size = 1, name = "classical4")
   classical5 = ClassicalRegister(size = 1, name = "classical5")
   classical6 = ClassicalRegister(size = 1, name = "classical6")
   classical7 = ClassicalRegister(size = 1, name = "classical7")
   classical8 = ClassicalRegister(size = 1, name = "classical8")
   classical9 = ClassicalRegister(size = 1, name = "classical9")
   classical10 = ClassicalRegister(size = 1, name = "classical10")
   classical11 = ClassicalRegister(size = 1, name = "classical11")
   classical12 = ClassicalRegister(size = 1, name = "classical12")
   
   qc = QuantumCircuit(env_qubit, code_qubits, ancillae_qubits, aux_qubits, bits, classical1, classical2, classical3, classical4, classical5, classical6, classical7, classical8, classical9, classical10, classical11, classical12) 
   
   qc.ry(theta = env_param, qubit = env_qubit)
   qc.barrier()
   
   qc.h(qubit = [code_qubits[0], code_qubits[3], code_qubits[6]])
   for i in range(0,9,3):
       
       qc.cx(control_qubit = code_qubits[i], target_qubit = [code_qubits[i+1], code_qubits[i+2]])
       
   qc.barrier()
   
   qc.cx(control_qubit = env_qubit, target_qubit = qubit_error_target)
   
   qc.barrier()
   
   if noisy_cnot is None:
       
       for i in range(0,2):
           
           qc.cx(control_qubit = code_qubits[i], target_qubit = ancillae_qubits[0])
   else:
       
       qc = qc.compose(noisy_cnot, qubits = [1,10])
       qc.cx(control_qubit = code_qubits[1], target_qubit = ancillae_qubits[0])
    
   for i in range(1,3):
       
       qc.cx(control_qubit = code_qubits[i], target_qubit = ancillae_qubits[1])
       
   qc.barrier()
   
   for i in range(3,5):
       
       qc.cx(control_qubit = code_qubits[i], target_qubit = ancillae_qubits[2])
   
   for i in range(4,6):
       
       qc.cx(control_qubit = code_qubits[i], target_qubit = ancillae_qubits[3])
       
   qc.barrier()
   
   for i in range(6,8):
       
       qc.cx(control_qubit = code_qubits[i], target_qubit = ancillae_qubits[4])
   
   for i in range(7,9):
       
       qc.cx(control_qubit = code_qubits[i], target_qubit = ancillae_qubits[5])
       
   qc.barrier()
   
   qc.mct(ancillae_qubits[:2], aux_qubits[0])
   qc.cx(control_qubit = aux_qubits[0], target_qubit = [ancillae_qubits[0], ancillae_qubits[1]])  
   qc.mct(ancillae_qubits[2:4], aux_qubits[1])
   qc.cx(control_qubit = aux_qubits[1], target_qubit = [ancillae_qubits[2], ancillae_qubits[3]])
   qc.mct(ancillae_qubits[4:6], aux_qubits[2])
   qc.cx(control_qubit = aux_qubits[2], target_qubit = [ancillae_qubits[4], ancillae_qubits[5]])
   
   qc.barrier()
   
   qc.h(qubit = code_qubits[:6])
   for i in range(0,6):
       
       qc.cx(control_qubit = code_qubits[i], target_qubit = ancillae_qubits[6])
   qc.h(qubit = code_qubits[:6])
   
   qc.barrier()
   
   qc.h(qubit = code_qubits[3:])
   for i in range(3,9):
       
       qc.cx(control_qubit = code_qubits[i], target_qubit = ancillae_qubits[7])
   qc.h(qubit = code_qubits[3:])
   
   qc.barrier()
   
   qc.mct(ancillae_qubits[6:], aux_qubits[3])
   qc.cx(control_qubit = aux_qubits[3], target_qubit = ancillae_qubits[6:])
   
   qc.barrier()
   
   qc.measure(qubit = ancillae_qubits[0], cbit = classical1)
   qc.measure(qubit = ancillae_qubits[1], cbit = classical2)
   qc.measure(qubit = ancillae_qubits[2], cbit = classical3)
   qc.measure(qubit = ancillae_qubits[3], cbit = classical4)
   qc.measure(qubit = ancillae_qubits[4], cbit = classical5)
   qc.measure(qubit = ancillae_qubits[5], cbit = classical6)
   qc.measure(qubit = ancillae_qubits[6], cbit = classical7)
   qc.measure(qubit = ancillae_qubits[7], cbit = classical8)
   qc.measure(qubit = aux_qubits[3], cbit = classical9)
   qc.measure(qubit = aux_qubits[0], cbit = classical10)
   qc.measure(qubit = aux_qubits[1], cbit = classical11)
   qc.measure(qubit = aux_qubits[2], cbit = classical12)
   
   qc.barrier()
   
   qc.x(qubit = code_qubits[0]).c_if(classical = classical1, val = 1)
   qc.x(qubit = code_qubits[1]).c_if(classical = classical10, val = 1)
   qc.x(qubit = code_qubits[2]).c_if(classical = classical2, val = 1)
   qc.x(qubit = code_qubits[3]).c_if(classical = classical3, val = 1)
   qc.x(qubit = code_qubits[4]).c_if(classical = classical11, val = 1)
   qc.x(qubit = code_qubits[5]).c_if(classical = classical4, val = 1)
   qc.x(qubit = code_qubits[6]).c_if(classical = classical5, val = 1)
   qc.x(qubit = code_qubits[7]).c_if(classical = classical12, val = 1)
   qc.x(qubit = code_qubits[8]).c_if(classical = classical6, val = 1)
   
   qc.barrier()
   
   qc.z(qubit = code_qubits[:3]).c_if(classical = classical7, val = 1)
   
   qc.barrier()
   
   qc.z(qubit = code_qubits[3:6]).c_if(classical = classical9, val = 1)
   
   qc.barrier()
   
   qc.z(qubit = code_qubits[6:]).c_if(classical = classical8, val = 1)
   
   qc.barrier()
   
   for i in range(0,9,3):
       
       qc.cx(control_qubit = code_qubits[i], target_qubit = [code_qubits[i+2], code_qubits[i+1]])        
   qc.h(qubit = [code_qubits[0], code_qubits[3], code_qubits[6]])
   
   qc.barrier()
   
   qc.measure(qubit = code_qubits, cbit = bits)
   
   return qc 