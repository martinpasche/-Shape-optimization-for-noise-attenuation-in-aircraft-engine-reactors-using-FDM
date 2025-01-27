import numpy as np
from math import *
import Mapping
import scipy

dx = 0.001
dy = 0.001

def BC_12(
    T, i, j, chi
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_12
    L = [
        T[i, j + 1]["flag"] == "fantome", # droite
        T[i, j - 1]["flag"] == "fantome", # gauche
        T[i + 1, j]["flag"] == "fantome", # bas 
        T[i - 1, j]["flag"] == "fantome", # haut
    ]
    l = sum(L)
    
    bord = [
            [
                (i, j),
                -complex(0,1)*chi*alphaR + chi*alphaI - 2*complex(0,1)*chi*alphaR*(M0 / (k0*dx))**2 + 2*chi*alphaI*(M0 / (k0*dx))**2 - 2*complex(0,1)*k0*M0,
                (-complex(0,1)*chi*alphaR + chi*alphaI + 2*complex(0,1)*chi*alphaR*(M0 / (k0*dx))**2 - 2*chi*alphaI*(M0 / (k0*dx))**2 - 2*complex(0,1)*k0*M0) * complex(0,1)

            ],
            [
                (i, j+1),
                (2*chi*alphaR*(M0 / k) + M0**2)/(2*dx) + (complex(0,1) * alphaR - alphaI) * chi*M0**2/(k0 * dx)**2 + M0**2 / (2*dx),
                ((2*chi*alphaR*(M0 / k) + M0**2)/(2*dx) - (complex(0,1) * alphaR - alphaI) * chi*M0**2/(k0 * dx)**2 + M0**2 / (2*dx)) * complex(0,1)

            ],
            [
                (i, j-1),
                -((2*chi*alphaR*(M0 / k) + M0**2)/(2*dx) + (complex(0,1) * alphaR - alphaI) * chi*M0**2/(k0 * dx)**2 + M0**2 / (2*dx)),
                (-((2*chi*alphaR*(M0 / k) + M0**2)/(2*dx) - (complex(0,1) * alphaR - alphaI) * chi*M0**2/(k0 * dx)**2 + M0**2 / (2*dx))) * complex(0,1)
            ],
            [
                (i+1, j),
                0,
                0
            ],
            [
                (i-1, j),
                0,
                0
            ]
    ] # donne l'équation des conditions sur le bord 1 ou 2
    
    D = ["in", "CGB", "CGH"]
            
    ret = [
        [
            (i, j),
            k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx)**2,
            (k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx)**2) * complex(0,1),
            T[i, j]["flag"] == "frontiere" and T[i, j] in D
        ],
        [
            (i, j+1),
            1 / dx**2 - (M0 / dx)**2 + complex(0,1)*k0*(M0 / dx),
            (1 / dx**2 - (M0 / dx)**2 + complex(0,1)*k0*(M0 / dx)) * complex(0,1),
            T[i, j]["flag"] == "frontiere" and T[i, j+1] in D
        ],
        [
            (i, j-1),
            1 / dx**2 - (M0 / dx)**2 - complex(0,1)*k0*(M0 / dx),
            (1 / dx**2 - (M0 / dx)**2 - complex(0,1)*k0*(M0 / dx)) * complex(0,1),
            T[i, j]["flag"] == "frontiere" and T[i, j-1] in D
        ],
        [
            (i+1, j),
            0,
            0,
            T[i, j]["flag"] == "frontiere" and T[i+1, j] in D
        ],
        [
            (i-1, j),
            0,
            0,
            T[i, j]["flag"] == "frontiere" and T[i-1, j] in D
        ]
    ] # équation sur omega, on remarque que [2] = i*[1]
   
    

    if l == 1:
        r = L.index(True)

        # ajout de la dérivée sur n
        
        if r < 2:
            bord[1][1] += - 1/2 * (-1)**r * 1/dx
            bord[1][2] += - 1/2 * (-1)**r * 1/dx
            
            bord[2][1] += 1/2 * (-1)**r * 1/dx
            bord[2][2] += 1/2 * (-1)**r * 1/dx
        else:
            bord[3][1] += 1/2 * (-1)**r * 1/dy
            bord[3][2] += 1/2 * (-1)**r * 1/dy
            
            bord[4][1] += - 1/2 * (-1)**r * 1/dy
            bord[4][2] += - 1/2 * (-1)**r * 1/dy
            
        
        # réinjection et suppression du point fantôme
        
        for i in range(5):
            if i != r+1:
                ret[i][1] += - bord[i][1] * (ret[r+1][1] / bord[r+1][1] + ret[r+1][2] / bord[r+1][2]) # à vérifier
                ret[i][2] += - bord[i][2] * (ret[r+1][1] / bord[r+1][1] + ret[r+1][2] / bord[r+1][2]) # à vérifier
        
        del ret[r + 1] # supprime le point fantôme de la liste ret
        return ret
    else:
        indices = [index for index, value in enumerate(L) if value == True]
        
        bord[1][1] += sqrt(2)/4 * (-1)**indices[0] * 1/dx
        bord[1][2] += sqrt(2)/4 * (-1)**indices[0] * 1/dx
        
        bord[2][1] += - sqrt(2)/4 * (-1)**indices[0] * 1/dx
        bord[2][2] += - sqrt(2)/4 * (-1)**indices[0] * 1/dx

        bord[3][1] += sqrt(2)/4 * (-1)**indices[1] * 1/dy
        bord[3][2] += sqrt(2)/4 * (-1)**indices[1] * 1/dy
        
        bord[4][1] += - sqrt(2)/4 * (-1)**indices[1] * 1/dy
        bord[4][2] += - sqrt(2)/4 * (-1)**indices[1] * 1/dy

            # réinjection et suppression du point fantôme
    
        for i in range(5):
            if i-1 not in indices:
                ret[i][1] += - bord[i][1] * (ret[indices[0]+1][1] / bord[indices[0]+1][1] + ret[indices[0]+1][2] / bord[indices[0]+1][2])
                ret[i][2] += - bord[i][2] * (ret[indices[0]+1][1] / bord[indices[0]+1][1] + ret[indices[0]+1][2] / bord[indices[0]+1][2])
                
                ret[i][1] += - bord[i][1] * (ret[indices[1]+1][1] / bord[indices[1]+1][1] + ret[indices[1]+1][2] / bord[indices[1]+1][2])
                ret[i][2] += - bord[i][2] * (ret[indices[1]+1][1] / bord[indices[1]+1][1] + ret[indices[1]+1][2] / bord[indices[1]+1][2])       
        del ret[indices[0] + 1] # supprime le point fantôme de la liste ret
        del ret[indices[1]] # supprime le point fantôme de la liste ret
        
        return ret
               
            
def BC_out(
    T, i, j, chi
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_12
    L = [
        T[i, j + 1]["flag"] == "fantome", # droite
        T[i, j - 1]["flag"] == "fantome", # gauche
        T[i + 1, j]["flag"] == "fantome", # bas
        T[i - 1, j]["flag"] == "fantome", # haut
    ]
    l = sum(L)
    
    bord = [
            [
                (i, j),
                complex(0,1)*k + 2*complex(0,1)*k*M0,
                -complex(0,1)*k + 2*complex(0,1)*k*M0
            ],
            [
                (i, j+1),
                -M0**2 / dx,
                -M0**2 / dx
            ],
            [
                (i, j-1),
                M0**2 / dx,
                M0**2 / dx
            ],
            [
                (i+1, j),
                0,
                0
            ],
            [
                (i-1, j),
                0,
                0
            ]
    ] # donne l'équation des conditions sur le bord out
    
    D = ["in", "CGB", "CGH"]
  
    ret = [
        [
            (i, j),
            k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx)**2,
            (k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx)**2) * complex(0,1),
            T[i, j]["flag"] == "frontiere" and T[i, j] in D
        ],
        [
            (i, j+1),
            1 / dx**2 - (M0 / dx)**2 + complex(0,1)*k0*(M0 / dx),
            (1 / dx**2 - (M0 / dx)**2 + complex(0,1)*k0*(M0 / dx)) * complex(0,1),
            T[i, j]["flag"] == "frontiere" and T[i, j+1] in D
        ],
        [
            (i, j-1),
            1 / dx**2 - (M0 / dx)**2 - complex(0,1)*k0*(M0 / dx),
            (1 / dx**2 - (M0 / dx)**2 - complex(0,1)*k0*(M0 / dx)) * complex(0,1),
            T[i, j]["flag"] == "frontiere" and T[i, j-1] in D
        ],
        [
            (i+1, j),
            0,
            0,
            T[i, j]["flag"] == "frontiere" and T[i+1, j] in D
        ],
        [
            (i-1, j),
            0,
            0,
            T[i, j]["flag"] == "frontiere" and T[i-1, j] in D
        ]
    ] # équation sur omega, on remarque que [2] = i*[1]
    

    if l == 1:
        r = L.index(True)
        
        # ajout de la dérivée sur n
        
        if r < 2:
            bord[1][1] += 1/2 * (-1)**r * 1/dx
            bord[1][2] += 1/2 * (-1)**r * 1/dx
            
            bord[2][1] += - 1/2 * (-1)**r * 1/dx
            bord[2][2] += - 1/2 * (-1)**r * 1/dx
        else:
            bord[3][1] += - 1/2 * (-1)**r * 1/dy
            bord[3][1] += - 1/2 * (-1)**r * 1/dy
            
            bord[4][1] += 1/2 * (-1)**r * 1/dy
            bord[4][1] += 1/2 * (-1)**r * 1/dy
        
         
        # réinjection et suppression du point fantôme
        
        for i in range(5):
            if i != r+1:
                ret[i][1] += - bord[i][1] * (ret[r+1][1] / bord[r+1][1] + ret[r+1][2] / bord[r+1][2]) # à vérifier
                ret[i][2] += - bord[i][2] * (ret[r+1][1] / bord[r+1][1] + ret[r+1][2] / bord[r+1][2]) # à vérifier
        
        del ret[r + 1] # supprime le point fantôme de la liste ret
        return ret
    else:
        indices = [index for index, value in enumerate(L) if value == True]
        
        bord[1][1] += sqrt(2)/4 * (-1)**indices[0] * 1/dx
        bord[1][2] += sqrt(2)/4 * (-1)**indices[0] * 1/dx
        
        bord[2][1] += - sqrt(2)/4 * (-1)**indices[0] * 1/dx
        bord[2][2] += - sqrt(2)/4 * (-1)**indices[0] * 1/dx

        bord[3][1] += - sqrt(2)/4 * (-1)**indices[1] * 1/dy
        bord[3][1] += - sqrt(2)/4 * (-1)**indices[1] * 1/dy
        
        bord[4][1] += sqrt(2)/4 * (-1)**indices[1] * 1/dy
        bord[4][1] += sqrt(2)/4 * (-1)**indices[1] * 1/dy

            # réinjection et suppression du point fantôme
    
        for i in range(5):
            if i-1 not in indices:
                ret[i][1] += - bord[i][1] * (ret[indices[0]+1][1] / bord[indices[0]+1][1] + ret[indices[0]+1][2] / bord[indices[0]+1][2]) # à vérifier
                ret[i][2] += - bord[i][2] * (ret[indices[0]+1][1] / bord[indices[0]+1][1] + ret[indices[0]+1][2] / bord[indices[0]+1][2]) # à vérifier
                
                ret[i][1] += - bord[i][1] * (ret[indices[1]+1][1] / bord[indices[1]+1][1] + ret[indices[1]+1][2] / bord[indices[1]+1][2]) # à vérifier
                ret[i][2] += - bord[i][2] * (ret[indices[1]+1][1] / bord[indices[1]+1][1] + ret[indices[1]+1][2] / bord[indices[1]+1][2]) # à vérifier        
        del ret[indices[0] + 1] # supprime le point fantôme de la liste ret
        del ret[indices[1]] # supprime le point fantôme de la liste ret
        
        for el in ret:
            valeur = [el[0], el[1].real, el[1].imag, el[2].real, el[2].imag, el[-1]]
            el = valeur
        
        return ret


def omega(
    T, i, j, chi
):
    D = ["in", "CGB", "CGH"]

    return [
        [
            (i, j),
            (k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx)**2).real,
            (k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx)**2).imag,
            ((k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx)**2) * complex(0,1)).real,
            ((k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx)**2) * complex(0,1)).imag,
            T[i, j]["flag"] == "frontiere" and T[i, j] in D
        ],
        [
            (i, j+1),
            (1 / dx**2 - (M0 / dx)**2 + complex(0,1)*k0*(M0 / dx)).real,
            (1 / dx**2 - (M0 / dx)**2 + complex(0,1)*k0*(M0 / dx)).imag,
            ((1 / dx**2 - (M0 / dx)**2 + complex(0,1)*k0*(M0 / dx)) * complex(0,1)).real,
            ((1 / dx**2 - (M0 / dx)**2 + complex(0,1)*k0*(M0 / dx)) * complex(0,1)).imag,
            T[i, j]["flag"] == "frontiere" and T[i, j+1] in D
        ],
        [
            (i, j-1),
            (1 / dx**2 - (M0 / dx)**2 - complex(0,1)*k0*(M0 / dx)).real,
            (1 / dx**2 - (M0 / dx)**2 - complex(0,1)*k0*(M0 / dx)).imag,
            ((1 / dx**2 - (M0 / dx)**2 - complex(0,1)*k0*(M0 / dx)) * complex(0,1)).real,
            ((1 / dx**2 - (M0 / dx)**2 - complex(0,1)*k0*(M0 / dx)) * complex(0,1)).imag,
            T[i, j]["flag"] == "frontiere" and T[i, j-1] in D
        ],
        [
            (i+1, j),
            0,
            0,
            0,
            0,
            T[i, j]["flag"] == "frontiere" and T[i+1, j] in D
        ],
        [
            (i-1, j),
            0,
            0,
            0,
            0,
            T[i, j]["flag"] == "frontiere" and T[i-1, j] in D
        ]
    ] 
    


def BC_in(
    T, i, j
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_in
    return [[(i, j), 1, 0, 1, 0, False]]


def BC_coin(
    T, i, j, chi
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_12
    L = [
        T[i, j + 1]["flag"] == "fantome", # droite
        T[i, j - 1]["flag"] == "fantome", # gauche
        T[i + 1, j]["flag"] == "fantome", # bas
        T[i - 1, j]["flag"] == "fantome", # haut
    ]
      
    indices = [index for index, value in enumerate(L) if value == True]
    
    if 1 in indices: # cas où on a un coin à sur le bord in à gauche 
        return [[(i, j), 1, 1, False]]
    
    
    # sinon on est dans le cas d'un coin entre 1,2 et out donc haut ou bas à droite
    
    bord_12 = [
            [
                (i, j),
                -complex(0,1)*chi*alphaR + chi*alphaI - 2*complex(0,1)*chi*alphaR*(M0 / (k0*dx))**2 + 2*chi*alphaI*(M0 / (k0*dx))**2 - 2*complex(0,1)*k0*M0,
                -complex(0,1)*chi*alphaR + chi*alphaI + 2*complex(0,1)*chi*alphaR*(M0 / (k0*dx))**2 - 2*chi*alphaI*(M0 / (k0*dx))**2 - 2*complex(0,1)*k0*M0
            ],
            [
                (i, j+1),
                (2*chi*alphaR*(M0 / k) + M0**2)/(2*dx) + (complex(0,1) * alphaR - alphaI) * chi*M0**2/(k0 * dx)**2 + M0**2 / (2*dx),
                (2*chi*alphaR*(M0 / k) + M0**2)/(2*dx) - (complex(0,1) * alphaR - alphaI) * chi*M0**2/(k0 * dx)**2 + M0**2 / (2*dx)

            ],
            [
                (i, j-1),
                -((2*chi*alphaR*(M0 / k) + M0**2)/(2*dx) + (complex(0,1) * alphaR - alphaI) * chi*M0**2/(k0 * dx)**2 + M0**2 / (2*dx)),
                -((2*chi*alphaR*(M0 / k) + M0**2)/(2*dx) - (complex(0,1) * alphaR - alphaI) * chi*M0**2/(k0 * dx)**2 + M0**2 / (2*dx))
            ],
            [
                (i+1, j),
                0,
                0
            ],
            [
                (i-1, j),
                0,
                0
            ]
    ] # donne l'équation des conditions sur le bord 1 ou 2
    
    
    bord_out = [
            [
                (i, j),
                complex(0,1)*k + 2*complex(0,1)*k*M0,
                -complex(0,1)*k + 2*complex(0,1)*k*M0
            ],
            [
                (i, j+1),
                -M0**2 / dx,
                -M0**2 / dx
            ],
            [
                (i, j-1),
                M0**2 / dx,
                M0**2 / dx
            ],
            [
                (i+1, j),
                0,
                0
            ],
            [
                (i-1, j),
                0,
                0
            ]
    ] # donne l'équation des conditions sur le bord out
   
    
    D = ["in", "CGB", "CGH"]

    ret = [
        [
            (i, j),
            k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx)**2,
            (k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx)**2) * complex(0,1),
            T[i, j]["flag"] == "frontiere" and T[i, j] in D
        ],
        [
            (i, j+1),
            1 / dx**2 - (M0 / dx)**2 + complex(0,1)*k0*(M0 / dx),
            (1 / dx**2 - (M0 / dx)**2 + complex(0,1)*k0*(M0 / dx)) * complex(0,1),
            T[i, j]["flag"] == "frontiere" and T[i, j+1] in D
        ],
        [
            (i, j-1),
            1 / dx**2 - (M0 / dx)**2 - complex(0,1)*k0*(M0 / dx),
            (1 / dx**2 - (M0 / dx)**2 - complex(0,1)*k0*(M0 / dx)) * complex(0,1),
            T[i, j]["flag"] == "frontiere" and T[i, j-1] in D
        ],
        [
            (i+1, j),
            0,
            0,
            T[i, j]["flag"] == "frontiere" and T[i+1, j] in D
        ],
        [
            (i-1, j),
            0,
            0,
            T[i, j]["flag"] == "frontiere" and T[i-1, j] in D
        ]
    ] # équation sur omega, on remarque que [2] = i*[1]
    
    
    # ajout de la dérivée sur n
    
    bord_out[1][1] += sqrt(2)/4 * 1/dx
    bord_out[1][2] += sqrt(2)/4 * 1/dx
    
    bord_out[2][1] += - sqrt(2)/4 * 1/dx
    bord_out[2][2] += - sqrt(2)/4 * 1/dx

    bord_12[3][1] += sqrt(2)/4 * (-1)**indices[1] * 1/dy
    bord_12[3][2] += sqrt(2)/4 * (-1)**indices[1] * 1/dy
    
    bord_12[4][1] += - sqrt(2)/4 * (-1)**indices[1] * 1/dy
    bord_12[4][2] += - sqrt(2)/4 * (-1)**indices[1] * 1/dy
    

    # réinjection et suppression du point fantôme
    
    for i in range(5):
        if i-1 not in indices:
            ret[i][1] += - bord_out[i][1] * (ret[indices[0]+1][1] / bord_out[indices[0]+1][1] + ret[indices[0]+1][2] / bord_out[indices[0]+1][2]) # à vérifier
            ret[i][2] += - bord_out[i][2] * (ret[indices[0]+1][1] / bord_out[indices[0]+1][1] + ret[indices[0]+1][2] / bord_out[indices[0]+1][2]) # à vérifier
             
            ret[i][1] += - bord_12[i][1] * (ret[indices[1]+1][1] / bord_12[indices[1]+1][1] + ret[indices[1]+1][2] / bord_12[indices[1]+1][2]) # à vérifier
            ret[i][2] += - bord_12[i][2] * (ret[indices[1]+1][1] / bord_12[indices[1]+1][1] + ret[indices[1]+1][2] / bord_12[indices[1]+1][2]) # à vérifier  
    
    del ret[indices[0] + 1] # supprime le point fantôme de la liste ret
    del ret[indices[1]] # supprime le point fantôme de la liste ret
    
    for el in ret:
        valeur = [el[0], el[1].real, el[1].imag, el[2].real, el[2].imag, el[-1]]
        el = valeur
        
    return ret


def vectorise_adj(
    i, j, p
):  # L'élement (i,j) a une partie réelle et une partie imaginaire, on les découpe en deux élements. La partie réelle est placée à la position 2pi + 2j et la partie imagainaire est placée à la position 2pi+2j+1
    return 2 * p * i + 2 * j


def matricise_adj(
    k, p
):  # matricise permet de retrouver les coordonnées du point correspondant à l'élément k du vecteur des valeurs dans la matrice de l'image. Appliquer Matricise à k_re et k_im renvoit le meme couple (i,j) qui correspont au meme point de l'image de base.
    return (k // (2 * p), (k % (2 * p)) // 2)


def matrice_to_csr(
    T, chi, b, z=1
):  # prend en entrée la matrice réelle avec les flags et le second membre
    n, p = T.shape[0], T.shape[1]

    # Initialisation des listes CSR
    A = []
    IA = []
    JA = []
    B = []  # second membre ajusté à cause des conditions de type dirichlet
    IB = (
        []
    )  # indice (i,j) éuivalent aux indices de B, est utile pour la fonction solution

    # Le premier élément de IA est toujours 0 pour la première rangée
    IA.append(0)

    for k in range(0, n * p):
        (i, j) = matricise_adj(k, p)
        if T[i, j] is not None and T[i, j]["flag"] == "EDP":
            L = omega(T, i, j, z)
            K = [
                index for index, sublist in enumerate(L) if sublist[3] == True
            ]  # donne les indices des voisins de i,j qui sont dans Gamma_in
            a = 0  # a est un compteur pour savoir si B[indice équivalent à la ligne i,j] existe déjà
            for r in range(
                0, len(L)
            ):  # L[r][0] parcourt les indices des voisins de i,j à prendre en compte
                (u, v) = L[r][0]
                q = vectorise_adj(u, v, p)
                if r in K:  # on voit si le voisin (u,v) est dans Gamma_in
                    if a == 0:  # B[equivalent i,j] n'existe pas encore
                        B.append(
                            b[k] - L[r][1] * b[q]
                        )  # on fait le relevé de dirichlet
                    else:
                        B[-1] += (
                            -L[r][1] * b[q]
                        )  # B[...] existe déjà et on lui ajoute la partie correspondante à un autre voisin (u,v)
                    a += 1
                else:  # si (u,v) n'est pas dans Gamma_in on ajoute sa valeur correspondante dans la matrice
                    A.append(L[r][1])
                    JA.append(q)
            if (
                a == 0
            ):  # si aucun voisin n'est dans Gamma_in B[...] on lui donne la valeur de b[k]
                B.append(b[k])
            IB.append((i, j))

        if T[i, j]["flag"] == "frontiere":
            if T[i, j]["BC"] == "in":
                L = BC_in(T, i, j, z)
            if T[i, j]["BC"] == "out":
                L = BC_out(T, i, j, z)
            if T[i, j]["BC"] == "12":
                L = BC_12(T, i, j, chi(i, j), z)
            if T[i, j]["BC"] == "Neumann":
                L = BC_12(T, i, j, chi(i, j), z)
            if T[i, j]["BC"] in ["CDB", "CDH", "CGH", "CGB"]:
                L = BC_coin(T, i, j, chi(i, j), z)
            K = [index for index, sublist in enumerate(L) if sublist[2] == True]
            a = 0
            for r in range(0, len(L)):
                (u, v) = L[r][0]
                q = vectorise_adj(u, v, p)
                if r in K:
                    if a == 0:  # B[equivalent i,j] n'existe pas encore
                        B.append(
                            b[k] - L[r][1] * b[q]
                        )  # on fait le relèvement de dirichlet
                    else:
                        B[-1] += (
                            -L[r][1] * b[q]
                        )  # B[...] existe déjà et on lui ajoute la partie correspondante à un autre voisin (u,v)
                    a += 1
                else:
                    A.append(L[r][1])
                    JA.append(q)
            if (
                a == 0
            ):  # si aucun voisin n'est dans Gamma_in B[...] on lui donne la valeur de b[k]
                B.append(b[k])
            IB.append((i, j))

        IA.append(len(A))  # Met à jour l'indicateur de rangée
    # Conversion des listes en tableaux NumPy
    A = np.array(A, dtype=complex)
    IA = np.array(IA, dtype=int)
    JA = np.array(JA, dtype=int)

    # Crée une matrice CSR
    csrmatrix = scipy.sparse.csr_matrix((A, JA, IA), shape=(n * p, n * p))

    # Supprimer les lignes et les colonnes nulles
    non_zero_rows = csrmatrix.getnnz(axis=1) > 0
    non_zero_cols = csrmatrix.getnnz(axis=0) > 0
    csrmatrix = csrmatrix[non_zero_rows][:, non_zero_cols]

    matrice_dense = csrmatrix.todense()

    return (csrmatrix, B, IB)


def matrice_to_csr(
    T, chi, b, z=1
):  # prend en entrée la matrice réelle avec les flags et le second membre
    n, p = T.shape[0], T.shape[1]

    # Initialisation des listes CSR
    A = []  # vecteur qui conserve les valeurs non nulles de la matrice creuse
    IA = []
    JA = (
        []
    )  # vecteur qui conserve les indices des colonnes des valeurs non nulles de la matrice
    B = []  # second membre ajusté à cause des conditions de type dirichlet
    IB = (
        []
    )  # indice (i,j) équivalent aux indices de B, est utile pour la fonction solution

    # Le premier élément de IA est toujours 0 pour la première rangée
    IA.append(0)
    k = 0
    while k < 2 * n * p:
        (i, j) = matricise_adj(k, p)
        if T[i, j] is not None and T[i, j]["flag"] == "EDP":
            L = omega(T, i, j, z)
            K = [
                index for index, sublist in enumerate(L) if sublist[3] == True
            ]  # donne les indices des voisins de i,j qui sont dans Gamma_in
            a = 0  # a est un compteur pour savoir si B[indice équivalent à la ligne i,j] existe déjà
            for r in range(
                0, len(L)
            ):  # L[r][0] parcourt les indices des voisins de i,j à prendre en compte
                (u, v) = L[r][0]
                q = vectorise_adj(u, v, p)
                if r in K:  # on voit si le voisin (u,v) est dans Gamma_in
                    if a == 0:  # B[equivalent i,j] n'existe pas encore
                        B.append(b[k] - L[r][1] * b[q])
                        B.append(
                            b[k + 1] - L[r][2] * b[q]
                        )  # on fait le relevé de dirichlet
                    else:  # PAS SUR LA ATTENTION PAS FORCEMENT JUSTE
                        B[-1] += -L[r][2] * b[q]
                        B[-2] += -L[r][1] * b[q]
                        # B[...] existe déjà et on lui ajoute la partie correspondante à un autre voisin (u,v)
                    a += 1
                else:  # si (u,v) n'est pas dans Gamma_in on ajoute sa valeur correspondante dans la matrice
                    A.append(L[r][1])
                    A.append(L[r][2])
                    JA.append(q)
                    JA.append(q + 1)
            if (
                a == 0
            ):  # si aucun voisin n'est dans Gamma_in B[...] on lui donne la valeur de b[k]
                B.append(b[k])
                B.append(b[k + 1])
            IB.append(
                (i, j)
            )  # PAS SUR A VOIR A MODIFIER PEUT ETRE DISTINGUER PARTIE REELLE ET IMAGINAIRE

        if T[i, j]["flag"] == "frontiere":
            if T[i, j]["BC"] == "in":
                L = BC_in(T, i, j)
            if T[i, j]["BC"] == "out":
                L = BC_out(T, i, j, z)
            if T[i, j]["BC"] == "12":
                L = BC_12(T, i, j, chi(i, j), z)
            if T[i, j]["BC"] == "Neumann":
                L = BC_12(T, i, j, chi(i, j), z)
            if T[i, j]["BC"] in ["CDB", "CDH", "CGH", "CGB"]:
                L = BC_coin(T, i, j, chi(i, j), z)
            K = [index for index, sublist in enumerate(L) if sublist[2] == True]
            a = 0
            for r in range(0, len(L)):
                (u, v) = L[r][0]
                q = vectorise_adj(u, v, p)
                if r in K:
                    if a == 0:  # B[equivalent i,j] n'existe pas encore
                        B.append(b[k] - L[r][1] * b[q])
                        B.append(
                            b[k + 1] - L[r][2] * b[q]
                        )  # on fait le relèvement de dirichlet
                    else:
                        B[-1] += -L[r][2] * b[q]
                        B[-2] += (
                            -L[r][1] * b[q]
                        )  # B[...] existe déjà et on lui ajoute la partie correspondante à un autre voisin (u,v)
                    a += 1
                else:
                    A.append(L[r][1])
                    A.append(L[r][2])
                    JA.append(q)
                    JA.append(q + 1)
            if (
                a == 0
            ):  # si aucun voisin n'est dans Gamma_in B[...] on lui donne la valeur de b[k]
                B.append(b[k])
                B.append(b[k + 1])
            IB.append(
                (i, j)
            )  # PAS SUR A VOIR A MODIFIER PEUT ETRE DISTINGUER PARTIE REELLE ET IMAGINAIRE

        IA.append(len(A))  # Met à jour l'indicateur de rangée
        k += 2  # A chaque itération on coupe la partie imaginaire et la partie réelle en 2 valeurs distinctes
    # Conversion des listes en tableaux NumPy
    A = np.array(A, dtype=complex)
    IA = np.array(IA, dtype=int)
    JA = np.array(JA, dtype=int)

    # Crée une matrice CSR
    csrmatrix = scipy.sparse.csr_matrix((A, JA, IA), shape=(n * p, 2 * n * p))

    # Supprimer les lignes et les colonnes nulles
    non_zero_rows = csrmatrix.getnnz(axis=1) > 0
    non_zero_cols = csrmatrix.getnnz(axis=0) > 0
    csrmatrix = csrmatrix[non_zero_rows][:, non_zero_cols]

    matrice_dense = csrmatrix.todense()

    return (csrmatrix, B, IB)


