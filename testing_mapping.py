import numpy as np
import fonction_analytique as anat
import matplotlib.pyplot as plt
from os import path
import Mapping as mapp

""" 
This file is used to debug the Mapping.reduction function.
The mapping.reduction function is used to define each node 
category in the domain, so the idea here is to visualize the
matrix that is returned by the reduction function to see if
the nodes are well defined.
"""

image_path = path.join("images", "Irregulier.png")

T = mapp.reduction(image_path)


""" T1 = np.array(list(map(lambda x: 65 if x == 2 else x, T.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 80 if x == 8 else x, T1.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 55 if x == 0 else x, T1.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 60 if x == 1 else x, T1.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 70 if x == 7 else x, T1.flatten()))).reshape(T.shape) """
T1 = np.array(list(map(lambda x: 80 if x == 72 else x, T.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 81 if x == 2 else x, T1.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 81 if x == 73 else x, T1.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 81 if x == 71 else x, T1.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 80 if x == 8 else x, T1.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 81 if x == 0 else x, T1.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 81 if x == 1 else x, T1.flatten()))).reshape(T.shape)
T1 = np.array(list(map(lambda x: 80 if x == 7 else x, T1.flatten()))).reshape(T.shape)

anat.VisualizeMatrix(T1)