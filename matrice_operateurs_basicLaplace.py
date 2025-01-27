import numpy as np
from parameters import *
import fonction_analytique as anat
from scipy.sparse import csr_matrix    
from typing import List, Union, Callable, Optional
from typing import Optional, Callable



""" 
This file contains the functions to build the matrices for the restrictions on space in 2 dimensions. 
However, this file is considering a SIMPLIFIED VERSION OF THE OPERATOR, which ONLY contains the Laplacian operator and
the identity operator. 
The restrictions are built following the matrix T, which is a matrix that contains the domain and the boundary conditions.
The matrix T is a matrix with the same shape as the domain, but with values that indicate the type of restriction in each cell.
The values are defined as follows:

0: Outside the domain
2: Inside the domain
7 or 70: Neumann. Left Boundary condition (wave going outside)
7_: Up and down boundary conditions. The number indicates the direction of the normal vector to the boundary
8_ or 8: Right Boundary condition (wave coming inside).
9_: Boundary Condition for Neumann and side BC simultaneously (For the corners)

1  2  3
4  X  5
6  7  8

Example:

73: Up right boundary condition for the upper side BC.
77: Down direction for the lower side BC.
"""


class CSRBuilder:
    
    """
    Class to build a csr matrix following a column-wise convention.
    
    add_value_XXX( row, col, value ): the value is added to a list of values
    with their corresponding row and column.
    
    build(): builds the csr matrix with the values added.
    """
    
    def __init__(self, T : np.ndarray = None, nodes_y : int = None, nodes_x : int = None):
        
        if T is not None:
            nodes_y, nodes_x = T.shape
        else:
            assert nodes_y is None or nodes_x is None
        
        self.nodes_y = nodes_y
        self.nodes_x = nodes_x
        self.data = []
        self.row_indexes = []
        self.col_indexes = []
        self.grid_size = nodes_y * nodes_x
        
    
    def get_node(self, row : int, col : int) -> int:
        return int(row + col * self.nodes_y)
        
    
    # this are the nodes of the big matrix that we want to build
    def add_value(self, row_big_matrix: int, col_big_matrix: int, value: complex):
        self.data.append(value)
        self.row_indexes.append(row_big_matrix)
        self.col_indexes.append(col_big_matrix)

    """ 
    The add value xxxx recieve the row and col of the node and the value to add
    Remember that the matrixes are built as column wise and upwards, e.i., 
    node - 1 = node down
    node + 1 = node up
    node + nodes_y = node right
    node - nodes_y = node left
    """
    
    def add_value_center(self, row: int, col: int, value: complex):
        node = self.get_node(row, col)
        node_value = self.get_node(row, col) 
        self.add_value(node, node_value, value)

    def add_value_right(self, row: int, col: int, value: complex):
        node = self.get_node(row, col)
        node_value = self.get_node(row, col + 1) 
        self.add_value(node, node_value, value)

    def add_value_left(self, row: int, col: int, value: complex):
        node = self.get_node(row, col)
        node_value = self.get_node(row, col - 1) 
        self.add_value(node, node_value, value)

    def add_value_up(self, row: int, col: int, value: complex):
        node = self.get_node(row, col)
        node_value = self.get_node(row + 1, col) 
        self.add_value(node, node_value, value)

    def add_value_down(self, row: int, col: int, value: complex):
        node = self.get_node(row, col)
        node_value = self.get_node(row - 1, col) 
        self.add_value(node, node_value, value)

    def build(self):
        data = np.array(self.data, dtype=complex)
        row_indexes = np.array(self.row_indexes)
        col_indexes = np.array(self.col_indexes)
        return csr_matrix((data, (row_indexes, col_indexes)), shape=(self.grid_size, self.grid_size))




def laplacien_csr(T : np.ndarray, metric : anat.Metric = None) -> csr_matrix:
    
    """ 
    Builds the laplacian restrictions for the domain.
    Adds the restriction only if inside the domain (value = 2 in T)
    
    Input:
    T: np.ndarray -> matrix with the domain
    metric: anat.Metric -> metric of the domain (sizes)
    
    Output:
    A: csr_matrix -> matrix with the restriction
    """
    
    assert T.ndim == 2  #matrix must be 2D
    csr_builder = CSRBuilder(T)
    
    print("Building laplacien csr")
    
    hx = metric.hx
    hy = metric.hy
        
    for i in range(csr_builder.nodes_y):
        for j in range(csr_builder.nodes_x):
            if T[i, j] == 2:
        
                csr_builder.add_value_center(i, j, -2.0 / (hx * hx) - 2.0 / (hy * hy))
                csr_builder.add_value_up(i, j, 1.0 / (hx * hx))
                csr_builder.add_value_down(i, j, 1.0 / (hx * hx))
                csr_builder.add_value_left(i, j, 1.0 / (hy * hy))
                csr_builder.add_value_right(i, j, 1.0 / (hy * hy))
                
                # The borders are dangerous, they might give errors if not treated properly
                if i == 0:
                    raise ValueError("Laplacian builder on border. row == 0. Possible error source from mapping or T.")
                if i == csr_builder.nodes_y - 1:
                    raise ValueError("Laplacian builder on border. row == csr_builder.nodes_y - 1. Possible error source from mapping or T.")
                if j == 0:
                    raise ValueError("Laplacian builder on border. column == 0. Possible error source from mapping or T.")
                if j == csr_builder.nodes_x - 1:
                    raise ValueError("Laplacian builder on border. column == csr_builder.nodes_x - 1. Possible error source from mapping or T.")
                    
    return csr_builder.build()


def identite_csr(T : np.ndarray, scalar: complex = None, metric : anat.Metric = None ) -> csr_matrix:
    
    """ 
    Builds the matrix for the restriction given by k * P , which
    will create an identity like matrix. The scalar is the value that
    will be multiplied by this matrix.
    Adds the restriction inside the domain (value = 2)
    
    Input:
    T: np.ndarray -> matrix with the domain
    scalar: complex -> scalar to multiply the matrix
    metric: anat.Metric -> metric of the domain (sizes)
    
    Output:
    A: csr_matrix -> matrix with the restriction
    """
    
    assert T.ndim == 2  #matrix must be 2D
    print("Building identite csr")
    
    if scalar is None:
        scalar = 1

    csr_builder = CSRBuilder(T)
    
    for i in range(csr_builder.nodes_y):
        for j in range(csr_builder.nodes_x):
            if T[i, j] == 2:
                csr_builder.add_value_center(i, j, scalar)
    return csr_builder.build()


def Dx2_csr (T: np.ndarray, scalar : complex = 1, metric : anat.Metric = None) -> csr_matrix:
    
    """ 
    Builds the matrix for the restriction given by -M^2 * d^2 P / dx^2
    Adds the restriction inside the domain (value = 2)
    
    Input:
    T: np.ndarray -> matrix with the domain
    scalar: complex -> scalar to multiply the matrix
    metric: anat.Metric -> metric of the domain (sizes)
    
    Output:
    A: csr_matrix -> matrix with the restriction
    """
    
    assert T.ndim == 2  #matrix must be 2D
    
    hx = metric.hx
    hy = metric.hy
    
    print("Building Dx2 csr")
    
    csr_builder = CSRBuilder(T)
    for i in range(csr_builder.nodes_y):
        for j in range(csr_builder.nodes_x):
            if T[i, j] == 2:
                csr_builder.add_value_center(i, j, -2.0 * scalar / hx ** 2)
                csr_builder.add_value_left(i, j, scalar / hx ** 2)
                csr_builder.add_value_right(i, j, scalar / hx ** 2)
                
                #The borders are dangerous, they might give errors if not treated properly
                
                if j == 0:
                    raise ValueError("Dx2 builder on border. column == 0. Possible error source from mapping or T.")
                if j == csr_builder.nodes_x - 1:
                    raise ValueError("Dx2 builder on border. column == csr_builder.nodes_x - 1. Possible error source from mapping or T.")
    return csr_builder.build()      


def Dx1_csr (T: np.ndarray, scalar : complex = 1, metric : anat.Metric = None) -> csr_matrix:
    
    """ 
    Builds the matrix A for the restriction given by -2 I k M d^1 P / dx
    
    Input:
    T: np.ndarray -> matrix with the domain
    scalar: complex -> scalar to multiply the matrix
    metric: anat.Metric -> metric of the domain (sizes)
    
    Output:
    A: csr_matrix -> matrix with the restriction
    """
    
    assert T.ndim == 2  #matrix must be 2D
    
    hx = metric.hx
    hy = metric.hy
    h = metric.h
    
    print("Building Dx1 csr")
    
    csr_builder = CSRBuilder(T)
    for i in range(csr_builder.nodes_y):
        for j in range(csr_builder.nodes_x):
            if T[i, j] == 2:
                csr_builder.add_value_left(i, j, -1 * scalar / (2 * hx))
                csr_builder.add_value_right(i, j, scalar / (2 * hx))
                
                #The borders are dangerous, they might give errors if not treated properly
                if j == 0:
                    raise ValueError("Dx1 builder on border. column == 0. Possible error source from mapping or T.")
                if j == csr_builder.nodes_x - 1:
                    raise ValueError("Dx1 builder on border. column == csr_builder.nodes_x - 1. Possible error source from mapping or T.")
            
    return csr_builder.build()  


def force_test(T: np.ndarray, f: Optional[Callable[[int, int], float]] = None, metric : anat.Metric = None) -> csr_matrix:
    
    """ 
    This adds the force to the restriction. Therefore, it is added 
    to the b vector in the equation Ax = b
    Adds the restriction inside the domain (value = 2) and in the Neumann BC (value = 7 or 70)
    
    input:
    T: np.ndarray -> matrix with the domain
    f: Optional[Callable[[int, int], float]] -> function that will be added to b
    metric: anat.Metric -> metric of the domain (sizes)
    
    output:
    b: csr_matrix -> vector with the force added
    """
    
    assert T.ndim == 2  #matrix must be 2D
    
    print("Building force test csr")
    
    nodes_y, nodes_x = T.shape
    b = np.zeros((nodes_y*nodes_x),dtype=complex)
    
    hx = metric.hx
    hy = metric.hy
    h = metric.h
    
    
    # The functions defines below should be defined in 
    # system_builder.py in SystemBuilder class with the method
    # set_functions(f, g)
    
    if f is None:
        f = lambda x, y: -2 * np.exp(complex(0, k0 * x))

    elif f is "Irregulier":
        Lx = metric.Lx
        Ly = metric.Ly
        
        f = lambda x,y:(1/(10000 * Lx**4 * y**4)) * np.exp(
        1j * k * x - (Ly * x) / (10 * Lx * y)
    ) * h * (
        6000 * Lx**3 * y**5 + Ly**5 * x * (x**2 + y**2) +
        200 * Ly * Lx**2 * y**4 * (x - 2j * k * y**2) -
        Ly**3 * y**2 * (x**3 - 400j * k * Lx**2 * y**2 + x * y**2) +
        10 * Ly**2 * Lx * y**3 * (-x**2 + 3 * y**2 + 2j * k * x * y**2) -
        10 * Ly**4 * Lx * y * (3 * x**2 + 3 * y**2 + 2j * k * x * y**2)
    )
    
    
    for i in range(nodes_y):
        for j in range(nodes_x):
            k = i + j * nodes_y
            x = j * hx
            y = i * hy
            
            if int(T[i,j]) == 2 or int(T[i,j]) in [7, 70]:
                b[k] += f(x, y)      
    return b


def BC_neumann_csr(T: np.ndarray, metric : anat.Metric = None) -> csr_matrix:
    
    """
    This add the restriction to the matrix A for the neumann boundary conditions (value = 70 in 
    matrix domain T). This corresponds to the left side of the domain.
    To apply this solution, we consider hx == hy
    
    Input:
    T: np.ndarray -> matrix with the domain
    metric: anat.Metric -> metric of the domain (sizes)
    
    Output:
    A: csr_matrix -> matrix with the restriction
    """
    
    print("Building BC neumann csr")
    
    hx = metric.hx
    hy = metric.hy
    h = metric.h
    
    csr_builder = CSRBuilder(T)
    
    for row in range(csr_builder.nodes_y):
        for col in range(csr_builder.nodes_x):
                
            if int(T[row,col]) in [7, 70]:
                csr_builder.add_value_center(row, col, complex( -4.0 / h ** 2 + k0 ** 2 , -2.0 * k0 / h))
                csr_builder.add_value_right(row, col, 2.0 / h ** 2)
                csr_builder.add_value_up(row, col, 1.0 / h ** 2)
                csr_builder.add_value_down(row, col, 1.0 / h ** 2)
                
    return csr_builder.build()


def BC_onde_csr(T, g : Optional[Callable[[int, int], float]] = None, metric : anat.Metric = None) -> Union[csr_matrix, csr_matrix]:
    
    """
    This add the restriction to the matrix A for the boundary conditions of the wave equation
    on the right hand side of the domain (value = 8 in matrix domain T). 
    Basically, this defines the incoming wave.
    That is why, we make each node get the value corresponding to the incoming wave by making
    1 the diagonal term of the matrix A and giving the value of the incoming wave to the vector b.
    
    The incoming wave is defined by the function P(x, Ly) = exp(i * k0 * x) * g(x, Ly)
    
    Input:
    T: np.ndarray -> matrix with the domain
    g: Optional[Callable[[int, int], float]] -> function that will be added to the matrix
    metric: anat.Metric -> metric of the domain (sizes)
    
    Output:
    A: csr_matrix -> matrix with the restriction
    b: csr_matrix -> vector b with the incoming wave added
    """
    
    print("Building BC onde csr")
    
    hx = metric.hx
    hy = metric.hy
    h = metric.h
      
    nodes_y, nodes_x = T.shape
    b = np.zeros((nodes_y * nodes_x), dtype=complex)
    csr_builder = CSRBuilder(T)
    
    if g is None:
        g = lambda x, y: 1
        
    for row in range(nodes_y):
        for col in range(nodes_x):
            k = row + col * nodes_y  # Indice de la cellule (i, j)
            x = col * hx
            y = row * hy

            if T[row, col] in list(range(81, 89)) or T[row, col] in [8, 8.0]:
                b[k] += np.exp(complex(0, k0 * x)) * g(x, y)
                csr_builder.add_value_center(row, col, 1)
        
    matrix = csr_builder.build()
    return matrix, b


def BC_dirichlet_csr(T, f : Optional[Callable[[int, int], float]] = None, metric : anat.Metric = None, bc_side = None) -> Union[csr_matrix, csr_matrix]:
    
    """ 
    This adds the restriction to the matrix A for the dirichlet boundary conditions
    (value = 71 - 79 in matrix domain T + 91 + 96). 
    This corresponds to the up and down side of the domain.
    To make the value 0 in the boundary, we add the value to the vector b and make
    1 the value in the diagonal of the matrix A to the correspoding node.
    
    Input:
    T: np.ndarray -> matrix with the domain
    (Unused)f: Optional[Callable[[int, int], float]] -> function that will be added to the matrix
    metric: anat.Metric -> metric of the domain (sizes)
    (Unused)bc_side: str -> side of the domain where the boundary condition is applied
    
    Output:
    A: csr_matrix -> matrix with the restriction
    b: csr_matrix -> vector b with 0 values in the boundary
    
    * The unsued variables are kept for compatibility with the other functions *
    """
    
    print("Building BC dirichlet csr")
    
    h = metric.h
    hx = metric.hx
    hy = metric.hy
    
    nodes_y, nodes_x = T.shape
    b = np.zeros((nodes_y * nodes_x), dtype=complex)
    csr_builder = CSRBuilder(T)
    
    for row in range(csr_builder.nodes_y):
        for col in range(csr_builder.nodes_x):
            k = row + col * nodes_y
            if T[row, col] in list(range(71, 79)) + [91, 96]:  
                
                b[k] += 0
                csr_builder.add_value_center(row, col, 1)         
                                
    matrix = csr_builder.build()        
    return matrix, b


def BC_up_down_csr(T, f : Optional[Callable[[int, int], float]] = None, bc_side :str = None ,metric : anat.Metric = None) -> Union[csr_matrix, csr_matrix]:
    
    """ 
    Adds the restriction to the matrix A for the BC for the upper and lower sides of the domain.
    (value = 71 - 79 in matrix domain T + 91 + 96). It can be dirichlet, neumann or robin.
    To apply this solution, we consider hx == hy.
    
    Input:
    T: np.ndarray -> matrix with the domain
    f: Optional[Callable[[int, int], float]] -> function that will be added to the vector b
    bc_side: str -> BC to be implemented
    metric: anat.Metric -> metric of the domain (sizes)
    
    Output:
    A: csr_matrix -> matrix with the restriction
    b: csr_matrix -> vector b with force values added
    """
    

    print("Building BC up and down sides csr")
    
    h = metric.h
    hx = metric.hx
    hy = metric.hy
    
    nodes_y, nodes_x = T.shape
    b = np.zeros((nodes_y * nodes_x), dtype=complex)
    csr_builder = CSRBuilder(T)
    
    if bc_side == "robin":
        eta = Z0 / Z
    elif bc_side == "neumann":
        eta = 0
     
    for row in range(csr_builder.nodes_y):
        for col in range(csr_builder.nodes_x):
            index = row + col * nodes_y
            x = col * hx
            y = row * hy
            
            """ Here im going to have the border conditions for dirichelt and neumann """
            
            if bc_side == "dirichlet":
                if T[row, col] in list(range(71, 79)) + [91, 96]:             
                    b[index] += 0
                    csr_builder.add_value_center(row, col, 1)
            
            
            elif bc_side == "neumann" or bc_side == "robin":
                if T[row, col] in list(range(71, 79)) + [91, 96]:  
                    b[index] += f(x, y)
                                       
                    
                if T[row, col] == 71 or T[row,col] == 91: #Upper left
                    csr_builder.add_value_center(row, col, -k0*complex(-4+h**2*k0**2 + 2*M0**2 + 8*M0*eta, -4*h*k0*eta) /h /complex(h*k0*(-1+2*M0*eta), 2*M0**2*eta) )
                    csr_builder.add_value_down(row, col, 2/h**2)
                    csr_builder.add_value_right(row, col, complex(2*h*k0*(-1+M0**2+2*M0*eta), -4*M0**2*eta ) / h**2 /complex(h*k0*(-1+2*M0*eta), 2*M0**2*eta )  )
                    
                elif T[row, col] == 72: # up
                    csr_builder.add_value_center(row, col, complex( (2*M0**2-4)/h**2 + k0**2, -4*eta*M0**2/(k0*h**3) - 2*k0*eta/h))
                    csr_builder.add_value_down(row, col, 2 / h**2)
                    csr_builder.add_value_right(row, col, complex( (1 - M0**2 - 2*eta*M0)/h**2, -k0*M0/h + 2*eta*M0**2/(k0 * h**3)))
                    csr_builder.add_value_left(row, col, complex( (1 - M0**2 + 2*eta*M0)/h**2, k0*M0/h + 2*eta*M0**2/(k0 * h**3)))
                    
                elif T[row, col] == 73: # Upper right
                    csr_builder.add_value_center(row, col, k0*complex(-4+h**2*k0**2 + 2*M0**2 - 8*M0*eta , -4*h*k0*eta ) / h / complex( h*(k0+2*k0*M0*eta) , -2*M0**2*eta ))
                    csr_builder.add_value_down(row, col, 2/h**2)
                    csr_builder.add_value_left(row, col, complex(-2*h*k0*(-1+M0**2-2*M0*eta), 4*M0**2*eta) / h**2 / complex( h*(k0 + 2*k0*M0*eta), -2*M0**2*eta ) )
                    
                elif T[row, col] == 74: #left
                    csr_builder.add_value_center(row, col, complex(h**3 * k0**3 + 2*h*k0*(M0**2 - 2 - 4*M0*eta), 2*h**2*k0**2*eta-4*M0**2*eta) / h**2 / complex(h*(k0+2*k0*M0*eta), 2*M0**2*eta))
                    csr_builder.add_value_up(row, col, 1 / h**2)
                    csr_builder.add_value_down(row, col, 1 / h**2)
                    csr_builder.add_value_right(row, col, -2*k0/h*(M0**2-1-2*M0*eta) / complex(h*(k0+2*k0*M0*eta), 2*M0**2*eta))
                    
                elif T[row, col] == 75: #right
                    csr_builder.add_value_center(row, col, complex(h**3*k0**3 + 2*h*k0*(-2 + M0**2 - 4*M0*eta, -2*h**2*k0**2*eta + 4*M0**2*eta)) / h**2 / complex(h* (k0 + 2*k0*M0*eta) , -2*M0**2*eta))
                    csr_builder.add_value_down(row, col, 1 / h**2)
                    csr_builder.add_value_up(row, col, 1/h**2)
                    csr_builder.add_value_left(row, col, -2*k0*(-1 + M0**2 - 2*M0*eta)/h/complex(h*(k0 + 2*k0*M0*eta) , -2*M0**2*eta ) )
                    
                elif T[row, col] == 76 or T[row, col] == 96: # Lower left
                    csr_builder.add_value_center(row, col, -k0 * complex(-4 + h**2*k0**2 + 2*M0**2 + 8*M0*eta, -4*h*k0*eta) / h / complex(h*k0*(-1 + 2*M0*eta), 2*M0**2*eta) )
                    csr_builder.add_value_up(row, col, 2 / h**2)
                    csr_builder.add_value_right(row, col, complex(2*h*k0*(-1 + M0**2 + 2*M0*eta), -4*M0**2*eta) / h**2 / complex(h*k0*(-1+2*M0*eta), 2*M0**2*eta))
                    
                elif T[row, col] == 77: # down
                    csr_builder.add_value_center(row, col, complex( (2*M0**2-4)/h**2 + k0**2, -4*eta*M0**2/(k0*h**3) - 2*k0*eta/h))
                    csr_builder.add_value_up(row, col, 2 / h**2)
                    csr_builder.add_value_right(row, col, complex( (1 - M0**2 - 2*eta*M0)/h**2, - k0*M0/h + 2*eta*M0**2/(k0 * h**3)))
                    csr_builder.add_value_left(row, col, complex( (1 - M0**2 + 2*eta*M0)/h**2, k0*M0/h + 2*eta*M0**2/(k0 * h**3)))
                
                elif T[row, col] == 78: # Lower right
                    csr_builder.add_value_center(row, col, k0*complex(-4 + h**2*k0**2 + 2*M0**2 - 8*M0*eta, -4*h*k0*eta ) / h / complex(h*(k0 + 2*k0*M0*eta), -2*M0**2*eta) )
                    csr_builder.add_value_up(row, col, 2 / h**2)
                    csr_builder.add_value_left(row, col, complex(-2*h*k0*(-1 + M0**2 + 2*M0*eta), 4*M0**2*eta) / h**2 / complex( h*(k0 + 2*k0*M0*eta), -2*M0**2*eta ))
                                                        
    matrix = csr_builder.build()        
    return matrix, b
