import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.providers.aer.noise import QuantumError
from typing import List

def controlled_rx_cnot(measurements: bool = False, noisy: QuantumError = None, params: List[float] = None) -> QuantumCircuit:
    """Creates a CNOT gate using a controlled rx gate.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.
        noisy (QuantumError, optional): A noisy that will added to the circuit. Defaults to None.
        params (List[float], optional): A list of coherent errors that will added to the circuit. Defaults to None.

    Raises:
        TypeError: If measurements is not a bool.
        TypeError: If noisy is not a QuantumError.
        TypeError: If noisy is not a QuantumError.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT gate.
    """
    if isinstance(measurements, bool):
        
        if measurements:
            
            qc = QuantumCircuit(2,2)
            if noisy is None:
                
                if params is None:
                    
                    qc.crx(theta = np.pi, control_qubit = 0, target_qubit = 1)
                else:
                    
                    qc.crx(theta = params[0], control_qubit = 0, target_qubit = 1)
                    
                qc.measure(qubit = [0,1], cbit = [0,1])
            else:
                
                if isinstance(noisy, QuantumError):
                    
                    qc.append(noisy, qargs = [0,1])
                    if params is None:
                        
                        qc.crx(theta = np.pi, control_qubit = 0, target_qubit = 1)
                    else:
                        
                        qc.crx(theta = params[0], control_qubit = 0, target_qubit = 1)
                        
                    qc.append(noisy, qargs = [0,1])
                    qc.measure(qubit = [0,1], cbit = [0,1])
                else:
                    
                    raise TypeError("The input is not a QuantumError!")
            
            return qc
        else:
            
            qc = QuantumCircuit(2)
            if noisy is None:
                
                if params is None:
                    
                    qc.crx(theta = np.pi, control_qubit = 0, target_qubit = 1)
                else:
                    
                    qc.crx(theta = params[0], control_qubit = 0, target_qubit = 1)
            else:
                
                if isinstance(noisy, QuantumError):
                    
                    qc.append(noisy, qargs = [0,1])
                    if params is None:
                    
                        qc.crx(theta = np.pi, control_qubit = 0, target_qubit = 1)
                    else:
                    
                        qc.crx(theta = params[0], control_qubit = 0, target_qubit = 1)
                    qc.append(noisy, qargs = [0,1])
                else:
                    
                    raise TypeError("The input is not a QuantumError!")
        
            return qc
    else:
        
        raise TypeError("The input is not a bool!")
    
def floating_gate_cnot(measurements: bool = False) -> QuantumCircuit:
    pass
    
def hadamard_cz_cnot(measurements: bool = False, noisy: QuantumError = None) -> QuantumCircuit:
    """Creates a CNOT gate using hadamard and CZ gates.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.
        noisy (QuantumError, optional): A noisy that will added to the circuit. Defaults to None.

    Raises:
        TypeError: If measurements is not a bool.
        TypeError: If noisy is not a QuantumError.
        TypeError: If noisy is not a QuantumError.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT gate.
    """
    if isinstance(measurements, bool):
        
        if measurements:
            
            qc = QuantumCircuit(2,2)
            if noisy is None:
                
                qc.h(qubit = 1)
                qc.cz(control_qubit = 0, target_qubit = 1)
                qc.h(qubit = 1)
                qc.measure(qubit = [0,1], cbit = [0,1])
            else:
                
                if isinstance(noisy, QuantumError):
                    
                    qc.append(noisy, qargs = [0,1])
                    qc.h(qubit = 1)
                    qc.cz(control_qubit = 0, target_qubit = 1)
                    qc.h(qubit = 1)
                    qc.append(noisy, qargs = [0,1])
                    qc.measure(qubit = [0,1], cbit = [0,1])
                else:
                    
                    raise TypeError("The input is not a QuantumError!")
            
            return qc
        else:
        
            qc = QuantumCircuit(2)
            if noisy is None:
                
                qc.h(qubit = 1)
                qc.cz(control_qubit = 0, target_qubit = 1)
                qc.h(qubit = 1)
            else:
                
                if isinstance(noisy, QuantumError):
                    
                    qc.append(noisy, qargs = [0,1])
                    qc.h(qubit = 1)
                    qc.cz(control_qubit = 0, target_qubit = 1)
                    qc.h(qubit = 1)
                    qc.append(noisy, qargs = [0,1])
                else:
                    
                    raise TypeError("The input is not a QuantumError!")
        
            return qc
    else:
        
        raise TypeError("The input is not a bool!")
    
def molmer_sorensen_cnot(measurements: bool = False, noisy: QuantumError = None, params: List[float] = None) -> QuantumCircuit:
    """Creates a CNOT gate using RY, RX and RXX gates.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.
        noisy (QuantumError, optional): A noisy that will added to the circuit. Defaults to None.
        params (List[float], optional): A list of coherent errors that will added to the circuit. Defaults to None.

    Raises:
        TypeError: If measurements is not a bool.
        TypeError: If noisy is not a QuantumError.
        TypeError: If noisy is not a QuantumError.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT gate.
    """
    if isinstance(measurements, bool):
        
        if measurements:
            
            qc = QuantumCircuit(2,2)
            if noisy is None:
                
                if params is None:
                    
                    qc.ry(theta = np.pi/2, qubit = 0)
                    qc.rxx(theta = np.pi/2, qubit1= 0, qubit2= 1)
                    qc.rx(theta = -np.pi/2, qubit = [0,1])
                    qc.ry(theta = -np.pi/2, qubit = 0)
                else:
                    
                    qc.ry(theta = params[0], qubit = 0)
                    qc.rxx(theta = params[1], qubit1= 0, qubit2= 1)
                    qc.rx(theta = params[2], qubit = [0,1])
                    qc.ry(theta = params[3], qubit = 0)
                    
                qc.measure(qubit = [0,1], cbit = [0,1])
            else:
                
                qc.append(noisy, qargs = [0,1])
                if params is None:
                    
                    qc.ry(theta = np.pi/2, qubit = 0)
                    qc.rxx(theta = np.pi/2, qubit1= 0, qubit2= 1)
                    qc.rx(theta = -np.pi/2, qubit = [0,1])
                    qc.ry(theta = -np.pi/2, qubit = 0)
                else:
                    
                    qc.ry(theta = params[0], qubit = 0)
                    qc.rxx(theta = params[1], qubit1= 0, qubit2= 1)
                    qc.rx(theta = params[2], qubit = [0,1])
                    qc.ry(theta = params[3], qubit = 0)
                    
                qc.append(noisy, qargs = [0,1])
                qc.measure(qubit = [0,1], cbit = [0,1])
            
            return qc
        else:
            
            qc = QuantumCircuit(2)
            if noisy is None:
                
                if params is None:
                    
                    qc.ry(theta = np.pi/2, qubit = 0)
                    qc.rxx(theta = np.pi/2, qubit1= 0, qubit2= 1)
                    qc.rx(theta = -np.pi/2, qubit = [0,1])
                    qc.ry(theta = -np.pi/2, qubit = 0)
                else:
                    
                    qc.ry(theta = params[0], qubit = 0)
                    qc.rxx(theta = params[1], qubit1= 0, qubit2= 1)
                    qc.rx(theta = params[2], qubit = [0,1])
                    qc.ry(theta = params[3], qubit = 0)
                    
            else:
                
                qc.append(noisy, qargs = [0,1])
                if params is None:
                    
                    qc.ry(theta = np.pi/2, qubit = 0)
                    qc.rxx(theta = np.pi/2, qubit1= 0, qubit2= 1)
                    qc.rx(theta = -np.pi/2, qubit = [0,1])
                    qc.ry(theta = -np.pi/2, qubit = 0)
                else:
                    
                    qc.ry(theta = params[0], qubit = 0)
                    qc.rxx(theta = params[1], qubit1= 0, qubit2= 1)
                    qc.rx(theta = params[2], qubit = [0,1])
                    qc.ry(theta = params[3], qubit = 0)
                    
                qc.append(noisy, qargs = [0,1]) 
        
            return qc
    else:
        
        raise TypeError("The input is not a bool!")
    
def standard_cnot(measurements: bool = False, noisy: QuantumError = None) -> QuantumCircuit:
    """Creates a quantum circuit with the standard CNOT gate.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.
        noisy (QuantumError, optional): A noisy that will added to the circuit. Defaults to False.

    Raises:
        TypeError: If measurements is not a bool.
        TypeError: If noisy is not a QuantumError.
        TypeError: If noisy is not a QuantumError.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT gate.
    """
    if isinstance(measurements, bool):
        
        if measurements:
            
            qc = QuantumCircuit(2,2)
            if noisy is None:
                
                qc.cx(control_qubit = 0, target_qubit = 1)
            else:
                
                if isinstance(noisy, QuantumError):
                    
                    qc.append(noisy, qargs = [0,1])
                    qc.cx(control_qubit = 0, target_qubit = 1)
                    qc.append(noisy, qargs = [0,1])
                else:
                    
                    raise TypeError("The input is not a QuantumError!")
            
            qc.measure(qubit = [0,1], cbit = [0,1])
            return qc
        else:
            
            qc = QuantumCircuit(2)
            if noisy is None:
                
                qc.cx(control_qubit = 0, target_qubit = 1)
            else:
                
                if isinstance(noisy, QuantumError):
                    
                    qc.append(noisy, qargs = [0,1])
                    qc.cx(control_qubit = 0, target_qubit = 1)
                    qc.append(noisy, qargs = [0,1])
                else:
                    
                    raise TypeError("The input is not a QuantumError!")
                
            return qc
    else:
        
        raise TypeError("The input is not a bool!")
    
def sqrt_swap_cnot(measurements: bool = False, noisy: QuantumError = None, params: List[float] = None) -> QuantumCircuit:
    """Creates a CNOT gate using the square root of the SWAP gate, RY, Z and RZ.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.
        noisy (QuantumError, optional): A noisy that will added to the circuit. Defaults to None.
        params (List[float], optional): A list of coherent errors that will added to the circuit. Defaults to None.

    Raises:
        TypeError: If measurements is not a bool.
        TypeError: If noisy is not a QuantumError.
        TypeError: If noisy is not a QuantumError.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT gate.
    """
    if isinstance(measurements, bool):
        
        sqrt_swap = np.array([[1,0,0,0],[0,0.5*(1+1j), 0.5*(1-1j),0],[0,0.5*(1-1j),0.5*(1+1j),0],[0,0,0,1]])
        if measurements:
            
            qc = QuantumCircuit(2,2)
            if noisy is None:
                
                if params is None:
                    
                    qc.ry(theta = np.pi/2, qubit = 1)
                    qc.unitary(sqrt_swap, qubits = [0,1])
                    qc.z(qubit = 0)
                    qc.unitary(sqrt_swap, qubits = [0,1])
                    qc.rz(phi = -np.pi/2, qubit = [0,1])
                    qc.ry(theta = -np.pi/2, qubit = 1)
                else:
                    
                    qc.ry(theta = params[0], qubit = 1)
                    qc.unitary(sqrt_swap, qubits = [0,1])
                    qc.z(qubit = 0)
                    qc.unitary(sqrt_swap, qubits = [0,1])
                    qc.rz(phi = params[1], qubit = [0,1])
                    qc.ry(theta = params[2], qubit = 1)
            else:
                
                if isinstance(noisy, QuantumError):
                    
                    qc.append(noisy, qargs = [0,1])
                    if params is None:
                    
                        qc.ry(theta = np.pi/2, qubit = 1)
                        qc.unitary(sqrt_swap, qubits = [0,1])
                        qc.z(qubit = 0)
                        qc.unitary(sqrt_swap, qubits = [0,1])
                        qc.rz(phi = -np.pi/2, qubit = [0,1])
                        qc.ry(theta = -np.pi/2, qubit = 1)
                    else:
                    
                        qc.ry(theta = params[0], qubit = 1)
                        qc.unitary(sqrt_swap, qubits = [0,1])
                        qc.z(qubit = 0)
                        qc.unitary(sqrt_swap, qubits = [0,1])
                        qc.rz(phi = params[1], qubit = [0,1])
                        qc.ry(theta = params[2], qubit = 1)
                    qc.append(noisy, qargs = [0,1])
                else:
                    
                    raise TypeError("The input is not a QuantumError!")
                
            qc.measure(qubit = [0,1], cbit = [0,1])
            return qc
        else:
            
            qc = QuantumCircuit(2)
            if noisy is None:
                
                if params is None:
                    
                    qc.ry(theta = np.pi/2, qubit = 1)
                    qc.unitary(sqrt_swap, qubits = [0,1])
                    qc.z(qubit = 0)
                    qc.unitary(sqrt_swap, qubits = [0,1])
                    qc.rz(phi = -np.pi/2, qubit = [0,1])
                    qc.ry(theta = -np.pi/2, qubit = 1)
                else:
                    
                    qc.ry(theta = params[0], qubit = 1)
                    qc.unitary(sqrt_swap, qubits = [0,1])
                    qc.z(qubit = 0)
                    qc.unitary(sqrt_swap, qubits = [0,1])
                    qc.rz(phi = params[1], qubit = [0,1])
                    qc.ry(theta = params[2], qubit = 1)
            else:
                
                if isinstance(noisy, QuantumError):
                    
                    qc.append(noisy, qargs = [0,1])
                    if params is None:
                    
                        qc.ry(theta = np.pi/2, qubit = 1)
                        qc.unitary(sqrt_swap, qubits = [0,1])
                        qc.z(qubit = 0)
                        qc.unitary(sqrt_swap, qubits = [0,1])
                        qc.rz(phi = -np.pi/2, qubit = [0,1])
                        qc.ry(theta = -np.pi/2, qubit = 1)
                    else:
                    
                        qc.ry(theta = params[0], qubit = 1)
                        qc.unitary(sqrt_swap, qubits = [0,1])
                        qc.z(qubit = 0)
                        qc.unitary(sqrt_swap, qubits = [0,1])
                        qc.rz(phi = params[1], qubit = [0,1])
                        qc.ry(theta = params[2], qubit = 1)
                    qc.append(noisy, qargs = [0,1])
                else:
                    
                    raise TypeError("The input is not a QuantumError!")
        
            return qc
    else:
        
        raise TypeError("The input is not a bool!")

    
def sqrt_xx_cnot(measurements: bool = False, noisy: QuantumError = None) -> QuantumCircuit:
    """Creates the circuit of square root of XX CNOT gate.

    Args:
        measurements (bool, optional): Put measurements gates in the circuit. Defaults to False.
        noisy (QuantumError, optional): A noisy that will added to the circuit. Defaults to None.

    Raises:
        TypeError: If measurements if not a bool.
        TypeError: If noisy is not a QuantumError. 
        TypeError: If noisy is not a QuantumError.

    Returns:
        QuantumCircuit: The circuit that represents a CNOT.
    """
    if isinstance(measurements, bool):
        
        identity_4_by_4 = np.eye(N=4)
        sigma_x = np.array([[0,1],[1,0]])
        sigma_xx = np.kron(sigma_x, sigma_x)
        sqrt_sigma_xx = ((1+1j)/np.sqrt(2))*((1/np.sqrt(2))*identity_4_by_4 + -1j*(1/np.sqrt(2))*sigma_xx)
        if measurements:
            
            qc = QuantumCircuit(2,2)
            if noisy is None:
                
                qc.x(qubit = 0)
                qc.h(qubit = 0)
                qc.unitary(sqrt_sigma_xx, qubits = [0,1])
                qc.h(qubit = 0)
                qc.sx(qubit = 1)
                qc.s(qubit = 0)
                qc.x(qubit = 0)
            else:
                
                if isinstance(noisy, QuantumError):
                    
                    qc.append(noisy, qargs = [0,1])
                    qc.x(qubit = 0)
                    qc.h(qubit = 0)
                    qc.unitary(sqrt_sigma_xx, qubits = [0,1])
                    qc.h(qubit = 0)
                    qc.sx(qubit = 1)
                    qc.s(qubit = 0)
                    qc.x(qubit = 0)
                    qc.append(noisy, qargs = [0,1])
                else:
                    
                    raise TypeError("The input is not a QuantumError!")
            
            qc.measure(qubit = [0,1], cbit = [0,1])
            return qc
        else:
            
            qc = QuantumCircuit(2)
            if noisy is None:
                
                qc.x(qubit = 0)
                qc.h(qubit = 0)
                qc.unitary(sqrt_sigma_xx, qubits = [0,1])
                qc.h(qubit = 0)
                qc.sx(qubit = 1)
                qc.s(qubit = 0)
                qc.x(qubit = 0)
            else:
                
                if isinstance(noisy, QuantumError):
                    
                    qc.append(noisy, qargs = [0,1])
                    qc.x(qubit = 0)
                    qc.h(qubit = 0)
                    qc.unitary(sqrt_sigma_xx, qubits = [0,1])
                    qc.h(qubit = 0)
                    qc.sx(qubit = 1)
                    qc.s(qubit = 0)
                    qc.x(qubit = 0)
                    qc.append(noisy, qargs = [0,1])
                else:
                    
                    raise TypeError("The input is not a QuantumError!")
            
            return qc
    else:
        
        raise TypeError("The input is not a bool!")