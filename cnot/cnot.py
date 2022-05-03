import numpy as np
from qiskit.circuit import ClassicalRegister, QuantumCircuit

def controlled_rx_cnot(measurements: bool = False) -> QuantumCircuit:
    """Creates a CNOT gate using a controlled rx gate.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.

    Raises:
        TypeError: If measurements is not a bool.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT gate.
    """
    if isinstance(measurements, bool):
        
        qc = QuantumCircuit(2)
        qc.crx(theta = np.pi, control_qubit = 0, target_qubit = 1)
        if measurements:
            
            bits = ClassicalRegister(2)
            qc.add_bits(bits = bits)
            qc.measure(qubit = [0,1], cbit = [0,1])
            
            return qc
        
        return qc
    else:
        
        raise TypeError("The input is not a bool!")
    
def floating_gate_cnot(measurements: bool = False) -> QuantumCircuit:
    pass
    
def hadamard_cz_cnot(measurements: bool = False) -> QuantumCircuit:
    """Creates a CNOT gate using hadamard and CZ gates.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.

    Raises:
        TypeError: If measurements is not a bool.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT gate.
    """
    if isinstance(measurements, bool):
        
        qc = QuantumCircuit(2)
        qc.h(qubit = 1)
        qc.cz(control_qubit = 0, target_qubit = 1)
        qc.h(qubit = 1)
        if measurements:
            
            bits = ClassicalRegister(2)
            qc.add_bits(bits = bits)
            qc.measure(qubit = [0,1], cbit = [0,1])
            
            return qc
        
        return qc
    else:
        
        raise TypeError("The input is not a bool!")
    
def molmer_sorensen_cnot(measurements: bool = False) -> QuantumCircuit:
    """Creates a CNOT gate using RY, RX and RXX gates.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.

    Raises:
        TypeError: If measurements is not a bool.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT gate.
    """
    if isinstance(measurements, bool):
        
        qc = QuantumCircuit(2)
        qc.ry(theta = np.pi/2, qubit = 0)
        qc.rxx(theta = np.pi/2, qubit1= 0, qubit2= 1)
        qc.rx(theta = -np.pi/2, qubit = [0,1])
        qc.ry(theta = -np.pi/2, qubit = 0)
        if measurements:
            
            bits = ClassicalRegister(2)
            qc.add_bits(bits = bits)
            qc.measure(qubit = [0,1], cbit = [0,1])
            
            return qc
        
        return qc
    else:
        
        raise TypeError("The input is not a bool!")
    
def standard_cnot(measurements: bool = False) -> QuantumCircuit:
    """Creates a quantum circuit with the standard CNOT gate.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.

    Raises:
        TypeError: If measurements is not a bool.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT gate.
    """
    if isinstance(measurements, bool):
        
        qc = QuantumCircuit(2)
        qc.cx(control_qubit = 0, target_qubit = 1)
        if measurements:
            
            bits = ClassicalRegister(2)
            qc.add_bits(bits = bits)
            qc.measure(qubit = [0,1], cbit = [0,1])
            
            return qc
        
        return qc
    else:
        
        raise TypeError("The input is not as a bool!")
    
def sqrt_swap_cnot(measurements: bool = False) -> QuantumCircuit:
    """Creates a CNOT gate using the square root of the SWAP gate, RY, Z and RZ.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.

    Raises:
        TypeError: If measurements is not a bool.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT gate.
    """
    if isinstance(measurements, bool):
        
        sqrt_swap = np.array([[1,0,0,0],[0,0.5*(1+1j), 0.5*(1-1j),0],[0,0.5*(1-1j),0.5*(1+1j),0],[0,0,0,1]])
        qc = QuantumCircuit(2)
        qc.ry(theta = np.pi/2, qubit = 1)
        qc.unitary(sqrt_swap, qubits = [0,1])
        qc.z(qubit = 0)
        qc.unitary(sqrt_swap, qubits = [0,1])
        qc.rz(phi = -np.pi/2, qubit = [0,1])
        qc.ry(theta = -np.pi/2, qubit = 1)
        if measurements:
            
            bits = ClassicalRegister(2)
            qc.add_bits(bits = bits)
            qc.measure(qubit = [0,1], bits = [0,1])
            
            return qc
        
        return qc
    else:
        
        raise TypeError("The input is not as a bool!")