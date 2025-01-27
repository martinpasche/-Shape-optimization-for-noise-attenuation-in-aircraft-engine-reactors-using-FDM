import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import inv
# import matrice_operateurs as operateur
from parameters import *
import matplotlib.pyplot as plt
import Mapping as mapp
import fonction_analytique as anat


# Matrice = np.array([[1,0,0,1,1,1,0,0],[1,1,1,1,0,1,1,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[1,1,1,1,0,1,1,0],[1,0,0,1,1,1,0,0],[0,0,0,0,0,0,0,0]])
Matrice = np.array([[0,0,0,0,0,0,0,1,1,1,1,0],
                    [1,1,1,1,0,0,0,1,0,0,1,0],
                    [1,0,0,1,1,1,1,1,0,0,1,0],
                    [1,0,0,0,0,0,0,0,0,0,1,0],
                    [1,0,0,0,0,0,0,0,0,0,1,0],
                    [1,0,0,0,0,0,0,0,0,0,1,0],
                    [1,0,0,0,0,0,0,0,0,0,1,0],
                    [1,1,1,1,1,1,1,1,0,0,1,0],
                    [0,0,0,0,0,0,0,1,1,1,1,0]])


Matrix = np.array([[0,0,0,0,0,0,0,1,1,1,1,0],
                    [1,1,1,1,0,0,0,1,0,0,1,0],
                    [1,0,0,1,1,1,1,1,0,0,1,0],
                    [1,0,0,0,0,0,0,0,0,0,1,0],
                    [1,1,1,1,1,1,1,1,1,0,1,0],
                    [1,0,0,0,0,0,0,0,1,0,1,0],
                    [1,1,1,0,1,1,1,0,1,0,1,0],
                    [1,0,1,0,1,0,1,0,1,0,1,0],
                    [1,0,1,1,1,0,1,1,1,0,1,0],
                    [1,0,0,0,0,0,0,0,0,0,1,0],
                    [1,1,1,1,1,1,1,0,0,0,1,0],
                    [0,0,0,0,0,0,1,1,1,1,1,0]])

# n, p = Matrice.shape
# j_min_u = 0
# j_min_d = 0
# j_max = p
# i_min_l = 1
# i_max_l = n-2
# i_min_r = i_min_l
# i_max_r = i_max_l

# borders = [[i_min_l,j_min_u]]
# upper_borders = [[i_max_l,j_min_u]]
# while j_min_u<j_max-2:
#     if Matrice[i_min_r,j_min_u+1]==1:
#         j_min_u+=1
#         borders.append([i_min_r,j_min_u])
#         Matrice[i_min_r,j_min_u]=5
#     elif Matrice[i_min_r+1,j_min_u]==1:
#         i_min_r+=1
#         borders.append([i_min_r,j_min_u])
#         Matrice[i_min_r,j_min_u]=5
#     elif Matrice[i_min_r-1,j_min_u]==1:
#         i_min_r-=1
#         borders.append([i_min_r,j_min_u])
#         Matrice[i_min_r,j_min_u]=5

# while j_min_d<j_max-2:

#     if Matrice[i_max_r,j_min_d+1]==1:
#         j_min_d+=1
#         upper_borders.append([i_max_r,j_min_d])
#         Matrice[i_max_r,j_min_d]=5
#     elif Matrice[i_max_r+1,j_min_d]==1:
#         i_max_r+=1
#         upper_borders.append([i_max_r,j_min_d])
#         Matrice[i_max_r,j_min_d]=5
#     elif Matrice[i_max_r-1,j_min_d]==1:
#         i_max_r-=1
#         upper_borders.append([i_max_r,j_min_d])
#         Matrice[i_max_r,j_min_d]=5
    

# # print("j u",j_min_u,'j d',j_min_d,"i min l",i_min_l,"i max l",i_max_l,"i min r",i_min_r,"i max r",i_max_r)
# # print(borders)
# # print(upper_borders)
# # print(Matrice)
# print(mapp.mapping_BC_opti_2(Matrix))

# image = mapp.color_to_flag_opti(mapp.mapping(mapp.png_to_rgb_matrix("réacteur2_mini.png")))

# for i in image :
#     for j in i:
#         if j!=0:
#             print("yeaaaah")
Matrix = np.zeros((100, 100))
print(Matrix.shape)
Lx = 100
Ly = 100
display = True

multiplier_4_increasing_resolution = 8
nodes_y, nodes_x = Matrix.shape
metric = anat.Metric(Lx, Ly, nodes_x, nodes_y)
T = anat.setting_matrix_domain(origin = "rect", nodes_x = Lx * multiplier_4_increasing_resolution + 1, nodes_y = Ly * multiplier_4_increasing_resolution + 1)
anat.plot_analytical_solution(Matrix,f = "irregulier",metric=metric,a=2)



