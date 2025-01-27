import numpy as np
import matplotlib.pyplot as plt
from parameters import *
import Mapping as mapp
from os import path
import os
import numpy as np
import pandas as pd
from typing import List, Union, Callable, Optional

""" 
In this file we can find all the complementary
functions that we can use in the other files.

In here we have defined the Metric class that
will be used to store the information about
the system.

We have also defined the setting_matrix_domain
function that will be used to set the domain
of the system.

We have also defined the plot_analytical_solution
function that will be used to plot the analytical
solution.

We have also defined the plot_numerical_real
function that will be used to plot the real part
of the numerical solution.

We have also defined the plot_numerical_img
function that will be used to plot the imaginary
part of the numerical solution.

We have also defined the erreurs function that
will be used to calculate the errors between the
analytical and numerical solutions. We are using norm
L2.

We have also defined the erreur_sur_domaine
function that will be used to calculate the difference
on the domain between the analytical and numerical solutions.

We have also defined the VisualizeSparseMatrix
function that will be used to visualize the sparse matrix.   <-- Very useful for debugging

We have also defined the VisualizeMatrix
function that will be used to visualize the matrix.  <-- Very useful for debugging

"""

class Metric:
    """
    Represents a metric for a given system.

    Attributes:
        Lx (float): The length of the system in the x-direction.
        Ly (float): The length of the system in the y-direction.
        nodes_x (int): The number of nodes in the x-direction.
        nodes_y (int): The number of nodes in the y-direction.
        dx (float): The spacing between nodes in the x-direction.
        dy (float): The spacing between nodes in the y-direction.
    """

    def __init__(self, Lx, Ly, nodes_x, nodes_y):
        self.Lx = Lx
        self.Ly = Ly
        self.nodes_x = nodes_x
        self.nodes_y = nodes_y
        self.dx = Lx / (nodes_x - 1)    # remember that we have nodes_x - 1 intervals
        self.dy = Ly / (nodes_y - 1)    # remember that we have nodes_y - 1 intervals

    @property
    def hx(self):
        return self.dx

    @property
    def hy(self):
        return self.dy

    @property
    def h(self):
        return self.dx

    def __str__(self):
        return "Lx = {}, Ly = {}, nodes_x = {}, nodes_y = {}, dx = {}, dy = {}".format(self.Lx, self.Ly, self.nodes_x, self.nodes_y, self.dx, self.dy)


def setting_matrix_domain(origin=None, nodes_x=None, nodes_y=None, image_path=None):
    """
    This function sets the matrix domain based on the provided origin shape.

    Parameters:
    - origin (str): The shape of the domain. Valid values are 
                        "square", "rectangle", "carre", "rectangulaire", 
                        "rect", "car", "squa", "rec", "rectangulaire", "squa", "rec" for a square or rectangle shape, 
                        "img" or "image" for an irregular shape based on an image, 
                        or "array" or "matrice" for a predefined array shape.
    - nodes_x (int): The number of nodes in the x-direction. Default is None.
    - nodes_y (int): The number of nodes in the y-direction. Default is None.
    - image_path (str): The path to the image file if origin is "img" or "image". Default is None.

    Returns:
    - T (numpy.ndarray): The matrix domain based on the provided origin shape.

    Raises:
    - Exception: If the provided origin shape is not valid or if image_path is not provided for "img" or "image" origin.

    """

    if origin is "square" or origin is "rectangle" or origin is "carre" or origin is "rectangulaire" \
        or "rect" in origin or "car" in origin or "squa" in origin or "rec" in origin \
        or "carre" in origin or "rectangle" in origin or "square" in origin \
        or "rectangulaire" in origin or "rect" in origin or "car" in origin \
        or "squa" in origin or "rec" in origin:
        
        if nodes_x is None:
            nodes_x = 60
        if nodes_y is None:
            nodes_y = 20       
            
        """ 
        In here we are building the matrix domain by hand.
        """
        
        T = np.ones((int(np.ceil(nodes_y)), int(np.ceil(nodes_x)))) * 2
        
        nodes_y, nodes_x = T.shape

        T[0, :] = 72
        T[nodes_y - 1, :] = 77
        T[ : , nodes_x- 1 ] = 8 #
        T[ : , 0] = 70 #

        T[0,0 ] = 91
        T[nodes_y - 1, 0] = 96
        
        T = np.flipud(T)
        
    
    elif "img" in origin or "image" in origin or "irregular" in origin:
        
        """ 
        In here we get the matrix domain from an image using the Mapping class.
        After building the matrix domain, we save it in a CSV file to avoid
        recomputing it every time.
        """
        
        if image_path is None:
            raise Exception("Please provide a valid path to the image")
        
        name_img = image_path.strip(".png").strip(".jpg").strip(".jpeg")
        
        try:
            os.mkdir("image_array")
        except FileExistsError:
                pass
        
        if not path.isfile( path.join("image_array", f"{name_img}.csv") ):
            T = mapp.reduction(path.join("images", image_path))
            T = np.flipud(T)
            np.savetxt(path.join("image_array", f"{name_img}.csv"), T, delimiter=",", fmt="%d")
        else:
            T = np.genfromtxt(path.join("image_array", f"{name_img}.csv"), delimiter=",")
        
        
        
    elif "array" in origin or "matrice" in origin:

        T = np.array([  [ 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
                        [74, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0],
                        [74, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0],
                        [74, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0],
                        [74, 2, 2, 2, 2, 2, 1, 1, 1, 1, 0],
                        [74, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0],
                        [74, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1],
                        [74, 2, 2, 2, 2, 2, 2, 2, 2, 2,85],
                        [74, 2, 2, 2, 2, 2, 2, 2, 2, 2,85],
                        [74, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1],
                        [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                        ])
    else: 
        raise Exception("Please provide a valid domain shape")
    return T


def plot_analytical_solution(u, a=1, f = None, metric = None):
    
    # this function can be cleaned up a bit. It is not clean code.
    
    """
    Plot the analytical solution of a given function.

    Parameters:
    - u: numpy.ndarray
        The input array representing the solution.
    - a: int, optional
        The plot type. Default is 1.
    - f: function, optional
        The function used to calculate the solution. Default is None.
    - metric: object, optional
        The metric object containing the dimensions of the plot. Default is None.

    Returns:
    - If a=3, returns Z_real and Z_imag as separate arrays.
    - If a=4 or a=5, returns the figure object and the axes objects.
    - Otherwise, returns the figure object and the axes objects.
    """
    
        
    nodes_y, nodes_x = u.shape
    
    Lx = metric.Lx
    Ly = metric.Ly
    
    
    # Normally, defining the function f is done 
    # in system_builder.py. However, we can define
    # it here if we want to plot the analytical
    # solution without solving the system.
    
    if f == None:
        f = lambda x, y: y * (Ly - y) * np.exp(complex(0, k * (x)))

    elif f == "Irregulier":
        f = lambda x,y:(y+Ly-1)/(10*(Lx-1)*(Ly-1))*\
            (y*y-(Ly-1)*y*(x/(10*(Lx-1))+1)+(Ly-1)*(Ly-1)*x/(10*(Lx-1)))\
            *np.exp(x*np.complex(-(Ly-1)/(10*(Lx-1))*(y+(Ly-1))/(y*y+(Ly-1)*y),k))\
            if y != 0 else 0

    Z = np.zeros((nodes_y, nodes_x), dtype=complex)
    for i in range(nodes_y):
        for j in range(nodes_x):
            x = j * metric.hx
            y = i * metric.hy
            
            Z[i,j]= f(x, y)


    # Obtenez la partie réelle et la partie imaginaire de Z
    Z_real = np.real(Z)
    Z_imag = np.imag(Z)
    

    # How to use np.linspace    
    """ 
    # Generate an array of evenly spaced values between start and stop
    x = np.linspace(start, stop, num)

    # Example usage
    start = 0
    stop = 10
    num = 5
    x = np.linspace(start, stop, num)
    print(x) """

    X, Y = np.meshgrid(np.linspace(0, metric.Lx, nodes_x), np.linspace(0, metric.Ly, nodes_y))

    if a==2:
        # Créez une figure 3D
        fig1 = plt.figure(figsize=(12, 5))

        # Graphique de la partie réelle
        ax1 = fig1.add_subplot(121, projection='3d')
        surf1 = ax1.plot_surface(X, Y, Z_real, cmap='coolwarm', edgecolor='none')
        ax1.set_title("Partie Réelle")
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax1.set_zlabel("Z_real")
        fig1.colorbar(surf1, ax=ax1)

        # Graphique de la partie imaginaire
        ax2 = fig1.add_subplot(122, projection='3d')
        surf2 = ax2.plot_surface(X, Y, Z_imag, cmap='coolwarm', edgecolor='none')
        ax2.set_title("Partie Imaginaire")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.set_zlabel("Z_imag")
        fig1.colorbar(surf2, ax=ax2)
        # Affichez la figure
        plt.show()
    
    elif a ==3:
        return Z_real,Z_imag
    
    elif a == 4:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Affichez la partie réelle avec une colormap coolwarm
        im1 = ax1.imshow(Z_real, cmap='coolwarm', origin='lower')
        ax1.set_title("Partie Réelle")
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        fig.colorbar(im1, ax=ax1)

        # Affichez la partie imaginaire avec une colormap coolwarm
        im2 = ax2.imshow(Z_imag, cmap='coolwarm', origin='lower')
        ax2.set_title("Partie Imaginaire")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        fig.colorbar(im2, ax=ax2)

        # Affichez la figure
        plt.show()
    
    elif a == 5:

        # Créez une figure 3D
        fig1 = plt.figure(figsize=(14, 8.5))

        # Graphique de la partie réelle
        ax1 = fig1.add_subplot(221, projection='3d')
        surf1 = ax1.plot_surface(X, Y, Z_real, cmap='coolwarm', edgecolor='none')
        ax1.set_title("Partie Réelle analytique")
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax1.set_zlabel("Pression")
        fig1.colorbar(surf1, ax=ax1)

        # Graphique de la partie réelle de u
        ax2 = fig1.add_subplot(223, projection='3d')
        surf2 = ax2.plot_surface(X, Y, u.real, cmap='coolwarm', edgecolor='none')
        ax2.set_title("Partie Réelle calculé")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.set_zlabel("Pression")
        fig1.colorbar(surf2, ax=ax2)

        # Graphique de la partie imaginaire
        ax2 = fig1.add_subplot(222, projection='3d')
        surf2 = ax2.plot_surface(X, Y, Z_imag, cmap='coolwarm', edgecolor='none')
        ax2.set_title("Partie Imaginaire analytique")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.set_zlabel("Pression")
        fig1.colorbar(surf2, ax=ax2)

        # Graphique de la partie réelle de u
        ax2 = fig1.add_subplot(224, projection='3d')
        surf2 = ax2.plot_surface(X, Y, u.imag, cmap='coolwarm', edgecolor='none')
        ax2.set_title("Partie imaginaire calculé")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.set_zlabel("Pression")
        fig1.colorbar(surf2, ax=ax2)

        return fig1,ax1,ax2


    else :
        # Créez une figure 3D
        fig1 = plt.figure(figsize=(12, 5))

        # Graphique de la partie réelle
        ax1 = fig1.add_subplot(121, projection='3d')
        surf1 = ax1.plot_surface(X, Y, Z_real, cmap='coolwarm', edgecolor='none')
        ax1.set_title("Partie Réelle")
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax1.set_zlabel("Z_real")
        fig1.colorbar(surf1, ax=ax1)

        # Graphique de la partie imaginaire
        ax2 = fig1.add_subplot(122, projection='3d')
        surf2 = ax2.plot_surface(X, Y, Z_imag, cmap='coolwarm', edgecolor='none')
        ax2.set_title("Partie Imaginaire")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.set_zlabel("Z_imag")
        fig1.colorbar(surf2, ax=ax2)
        # Ajouter des titres aux graphiques de la deuxième figure
        ax1.set_title("partie réelle analityque")
        ax2.set_title("partie imaginaire analityque")
        return fig1, ax1, ax2


def plot_numerical_real(u, fig, metric : Metric = None):
    """
    Plot the real part of a numerical solution in 3D.

    Parameters:
    u (ndarray): The numerical solution.
    fig (Figure): The figure object to plot on.
    metric (Metric, optional): The metric object containing the dimensions of the grid. Defaults to None.

    Returns:
    ax1 (Axes3D): The 3D axes object.
    fig (Figure): The updated figure object.
    """    
    
    # Définir les dimensions de votre grille 2D
    nodes_y, nodes_x = u.shape

    # Créer les grilles X et Y de manière appropriée
    X, Y = np.meshgrid(np.linspace(0, metric.Lx, nodes_x), np.linspace(0, metric.Ly, nodes_y))

    # Afficher la solution u en 3D (partie réelle)
    ax1 = fig.add_subplot(121, projection='3d')
    surf1 = ax1.plot_surface(X, Y, np.real(u), cmap='viridis', edgecolor='none')
    fig.colorbar(surf1, ax=ax1)
    
    #ax1.set_box_aspect([np.ptp(X), np.ptp(Y), min(np.ptp(X), np.ptp(Y))])
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('u_real')
    ax1.set_title("Partie Réelle calcué")
    
    return ax1, fig


def plot_numerical_img(u, fig, metric = None):
    """
    Plot the numerical image of the solution u in 3D (partie imaginaire).

    Parameters:
    u (numpy.ndarray): The solution array.
    fig (matplotlib.figure.Figure): The figure object to plot on.
    metric (object, optional): The metric object containing the dimensions of the grid. Defaults to None.

    Returns:
    matplotlib.axes._subplots.Axes3DSubplot: The 3D subplot object.
    matplotlib.figure.Figure: The updated figure object.
    """
    
    # Définir les dimensions de votre grille 2D
    nodes_y, nodes_x = u.shape

    # Créer les grilles X et Y de manière appropriée
    X, Y = np.meshgrid(np.linspace(0, metric.Lx, nodes_x), np.linspace(0, metric.Ly, nodes_y))

    # Afficher la solution u en 3D (partie imaginaire)
    ax2 = fig.add_subplot(122, projection='3d')
    surf2 = ax2.plot_surface(X, Y, np.imag(u), cmap='viridis', edgecolor='none')
    fig.colorbar(surf2, ax=ax2)
    
    #ax2.set_box_aspect([np.ptp(X), np.ptp(Y), min(np.ptp(X), np.ptp(Y))])
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('u_imag')
    ax2.set_title("Partie Imaginaire calculé")
    
    return ax2, fig



def erreurs(u, f = None, metric : Metric = None) -> Union[float, float]:
    """
    Calculate the errors between the analytical solution and the numerical solution.
    The errors are calculated using the L2 norm.

    Parameters:
    - u: numpy.ndarray
        The numerical solution.
    - f: function, optional
        The analytical solution function. Default is None.
    - metric: Metric, optional
        The metric object containing information about the grid. Default is None.

    Returns:
    - Real_errors: float
        The percentage error in the real part of the solution.
    - imag_errors: float
        The percentage error in the imaginary part of the solution.
    """
    
    Z = np.zeros((metric.nodes_y, metric.nodes_x), dtype=complex)
    for i in range(metric.nodes_y):
        for j in range(metric.nodes_x):
            x = j * metric.hx
            y = i * metric.hy
            
            Z[i,j]= f(x, y)

    # Obtenez la partie réelle et la partie imaginaire de Z
    Z_real = np.real(Z)
    Z_imag = np.imag(Z)
    
    nx,ny = u.shape
    
    hx = metric.hx
    hy = metric.hy
    
    Real_errors = 0
    imag_errors = 0
    real_anat = 0
    imag_anat = 0

    for j in range(ny):
        for i in range(nx):
            Real_errors += abs(Z_real[i,j]-np.real(u[i,j]))**2*hx*hy
            imag_errors +=abs(Z_imag[i,j]-np.imag(u[i,j]))**2*hx*hy
            real_anat += Z_real[i,j]**2*hx*hy
            imag_anat += Z_imag[i,j]**2*hx*hy
    
    epsilon = 1e-10
    if abs(real_anat) <= epsilon:
        Real_errors = -1
    else:
        Real_errors = np.sqrt(Real_errors)/np.sqrt(real_anat)*100
        
    if abs(imag_anat) <= epsilon:
        imag_errors = -1
    else:
        imag_errors = np.sqrt(imag_errors)/np.sqrt(imag_anat)*100
    
    return Real_errors,imag_errors



def erreur_sur_domaine(u : np.ndarray, erreur_real, erreur_img, f = None, metric : Metric = None) -> Union[np.ndarray, np.ndarray]:
    """
    Calculate the error on the domain between the given function u and the analytical function f.
    We say error, but it is the difference between the two functions.

    Parameters:
    u (complex array): The given function on the domain.
    (unused) erreur_real (float array): The array to store the error in the real part of the function.
    (unused) erreur_img (float array): The array to store the error in the imaginary part of the function.
    f (function, optional): The analytical function. Defaults to None.
    metric (Metric, optional): The metric object containing information about the domain. Defaults to None.

    Returns:
    tuple: A tuple containing the error in the real part and the error in the imaginary part of the function.
    """

    # Calculate the analytical function values on the domain
    Z = np.zeros((metric.nodes_y, metric.nodes_x), dtype=complex)
    for i in range(metric.nodes_y):
        for j in range(metric.nodes_x):
            x = j * metric.hx
            y = i * metric.hy
            
            Z[i,j]= f(x, y)

    # Obtain the real and imaginary parts of Z
    Z_real = np.real(Z)
    Z_imag = np.imag(Z)
    
    # Calculate the difference between the given function and the analytical function
    u_diff_real = abs(Z_real - np.real(u)) 
    u_diff_imag = abs(Z_imag - np.imag(u)) 

    return u_diff_real, u_diff_imag



def VisualizeSparseMatrix(sparse_matrix, show=True, title=None):
    """
    Visualizes a sparse matrix by plotting its real and imaginary parts.

    Parameters:
    sparse_matrix (scipy.sparse.spmatrix): The sparse matrix to visualize.
    show (bool, optional): Whether to display the plot. Defaults to True.
    title (str, optional): The title of the plot. Defaults to None.
    """

    # Extract real and imaginary parts
    real_Dx = np.real(sparse_matrix.toarray())
    imag_Dx = np.imag(sparse_matrix.toarray())

    # Create subplots
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    if title:
        fig.suptitle(title)

    # Plot real part
    cax0 = axs[0].imshow(real_Dx, cmap='viridis')
    axs[0].set_title('Real part')
    fig.colorbar(cax0, ax=axs[0])

    # Plot imaginary part
    cax1 = axs[1].imshow(imag_Dx, cmap='viridis')
    axs[1].set_title('Imaginary part')
    fig.colorbar(cax1, ax=axs[1])

    if show:
        plt.show()
    
    
def VisualizeMatrix(matrix, minmax=None):
    """
    Visualizes a complex matrix by plotting its real and imaginary parts.

    Parameters:
    matrix (numpy.ndarray): The complex matrix to visualize.
    minmax (tuple, optional): A tuple containing the minimum and maximum values for the color scale. 
                              If not provided, the color scale will be determined automatically.

    Returns:
    None
    """

    # Extract real and imaginary parts
    real_Dx = np.real(matrix)
    imag_Dx = np.imag(matrix)

    if minmax:
        minimum = minmax[0]
        maximum = minmax[1]

        # Create subplots
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

        # Plot real part
        cax0 = axs[0].imshow(real_Dx, cmap='viridis', vmin=minimum, vmax=maximum)
        axs[0].set_title('Real part')
        fig.colorbar(cax0, ax=axs[0])

        # Plot imaginary part
        cax1 = axs[1].imshow(imag_Dx, cmap='viridis', vmin=minimum, vmax=maximum)
        axs[1].set_title('Imaginary part')
        fig.colorbar(cax1, ax=axs[1])

        plt.show()

    else:
        # Create subplots
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

        # Plot real part
        cax0 = axs[0].imshow(real_Dx, cmap='viridis')
        axs[0].set_title('Real part')
        fig.colorbar(cax0, ax=axs[0])

        # Plot imaginary part
        cax1 = axs[1].imshow(imag_Dx, cmap='viridis')
        axs[1].set_title('Imaginary part')
        fig.colorbar(cax1, ax=axs[1])

        try:
            # Save the real part of the matrix as an Excel file
            excel_file = path.join("csv", "visualize_matrix.xlsx")
            df_real = pd.DataFrame(np.real(real_Dx))
            df_real.to_excel(excel_file, index=False)

            # Save the real part of the matrix as a CSV file
            csv_file = path.join("csv", "visualize_matrix.csv")
            df_real.to_csv(csv_file, index=False)
        except Exception as e:
            print(f"Error while saving the matrix as an Excel file: {e}")

        plt.show()
        
        
        


        