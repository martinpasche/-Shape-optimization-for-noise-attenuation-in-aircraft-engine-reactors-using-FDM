from scipy.sparse.linalg import spsolve
import numpy as np
import Mapping
import importlib
from math import *
import scipy.sparse
import cmath
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


def f_dx(k):
    return 0.01 / k


def f_dy(k):
    return 0.01 / k


k = 1
k0 = 1
M0 = 0.3
Z0 = 0.0004
Z = complex(1.1, -0.1)
alpha = k0 * Z0 / Z
alphaR = alpha.real
alphaI = alpha.imag


def BC_12(
    T, i, j, chi, z=1
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_12
    dx = f_dx(z)
    dy = f_dy(z)
    L = [
        T[i, j + 1]["flag"] == "fantome",  # droite
        T[i, j - 1]["flag"] == "fantome",  # gauche
        T[i + 1, j]["flag"] == "fantome",  # bas
        T[i - 1, j]["flag"] == "fantome",  # haut
    ]
    l = sum(L)

    bord = [
        [
            (i, j),
            -complex(0, 1) * chi * alphaR
            + chi * alphaI
            - 2 * complex(0, 1) * chi * alphaR * (M0 / (k0 * dx)) ** 2
            + 2 * chi * alphaI * (M0 / (k0 * dx)) ** 2
            - 2 * complex(0, 1) * k0 * M0,
            (-complex(0, 1) * chi * alphaR
            + chi * alphaI
            + 2 * complex(0, 1) * chi * alphaR * (M0 / (k0 * dx)) ** 2
            - 2 * chi * alphaI * (M0 / (k0 * dx)) ** 2
            - 2 * complex(0, 1) * k0 * M0) * complex(0, 1),
        ],
        [
            (i, j + 1),
            (2 * chi * alphaR * (M0 / k0) + M0**2+2*chi*alphaI*M0/k0) / (2 * dx)
            + (complex(0, 1) * alphaR - alphaI) * chi * M0**2 / (k0 * dx) ** 2
            + M0**2 / (2 * dx),
            ((2 * chi * alphaR * (M0 / k0) + M0**2+2*chi*alphaI*M0/k0) / (2 * dx)
            - (complex(0, 1) * alphaR - alphaI) * chi * M0**2 / (k0 * dx) ** 2
            + M0**2 / (2 * dx)) * complex(0, 1),
        ],
        [
            (i, j - 1),
            -(
                (2 * chi * alphaR * (M0 / k0) + M0**2+2*chi*alphaI*M0/k0) / (2 * dx)
                + (complex(0, 1) * alphaR - alphaI) * chi * M0**2 / (k0 * dx) ** 2
                + M0**2 / (2 * dx)
            ),
            -(
                (2 * chi * alphaR * (M0 / k0) + M0**2+2*chi*alphaI*M0/k0) / (2 * dx)
                - (complex(0, 1) * alphaR - alphaI) * chi * M0**2 / (k0 * dx) ** 2
                + M0**2 / (2 * dx)
            ) * complex(0, 1),
        ],
        [(i + 1, j), 0, 0],
        [(i - 1, j), 0, 0],
    ]  # donne l'équation des conditions sur le bord 1 ou 2
    D = ["in", "CGB", "CGH"]

    ret = [
        [
            (i, j),
            k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2,
            (k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2) * complex(0, 1),
            T[i, j]["flag"] == "frontiere" and T[i, j] in D,
        ],
        [
            (i, j + 1),
            1 / dx**2 - (M0 / dx) ** 2 + complex(0, 1) * k0 * (M0 / dx),
            (1 / dx**2 - (M0 / dx) ** 2 + complex(0, 1) * k0 * (M0 / dx))
            * complex(0, 1),
            T[i, j]["flag"] == "frontiere" and T[i, j + 1] in D,
        ],
        [
            (i, j - 1),
            1 / dx**2 - (M0 / dx) ** 2 - complex(0, 1) * k0 * (M0 / dx),
            (1 / dx**2 - (M0 / dx) ** 2 - complex(0, 1) * k0 * (M0 / dx))
            * complex(0, 1),
            T[i, j]["flag"] == "frontiere" and T[i, j - 1] in D,
        ],
        [(i + 1, j), 0, 0, T[i, j]["flag"] == "frontiere" and T[i + 1, j] in D],
        [(i - 1, j), 0, 0, T[i, j]["flag"] == "frontiere" and T[i - 1, j] in D],
    ]  # équation sur omega, on remarque que [2] = i*[1]

    if l == 1:
        r = L.index(True)

        # ajout de la dérivée sur n

        if r < 2:
            bord[1][1] += -1 / 2 * (-1) ** r * 1 / dx
            bord[1][2] += (-1 / 2 * (-1) ** r * 1 / dx)*complex(0,1)

            bord[2][1] += 1 / 2 * (-1) ** r * 1 / dx
            bord[2][2] += (1 / 2 * (-1) ** r * 1 / dx)*complex(0,1)
        else:
            bord[3][1] += 1 / 2 * (-1) ** r * 1 / dy
            bord[3][2] += (1 / 2 * (-1) ** r * 1 / dy)*complex(0,1)

            bord[4][1] += -1 / 2 * (-1) ** r * 1 / dy
            bord[4][2] += (-1 / 2 * (-1) ** r * 1 / dy)*complex(0,1)

        # réinjection et suppression du point fantôme

        for i in range(5):
            if i != r + 1:
                A=np.array([[bord[r+1][1].real,bord[r+1][2].real],[bord[r+1][1].imag,bord[r+1][2].imag]])
                B=-np.array([[bord[i][1].real,bord[i][2].real],[bord[i][1].imag,bord[i][2].imag]])
                C=np.linalg.inv(A)@B
                ret[i][1]+=ret[r+1][1]*C[0,0]+ret[r+1][2]*C[1,0]
                ret[i][2]+=ret[r+1][1]*C[0,1]+ret[r+1][2]*C[1,1]
                # ret[i][1] += -bord[i][1] * (
                #     ret[r + 1][1] / bord[r + 1][1] + ret[r + 1][2] / bord[r + 1][2]
                # )  # à vérifier
                # ret[i][2] += -bord[i][2] * (
                #     ret[r + 1][1] / bord[r + 1][1] + ret[r + 1][2] / bord[r + 1][2]
                # )  # à vérifier

        del ret[r + 1]  # supprime le point fantôme de la liste ret
        ret2 = ret.copy()
        for indice in range(len(ret)):
            el = ret[indice]
            valeur = [
                el[0],
                complex(el[1]).real,
                complex(el[1]).imag,
                complex(el[2]).real,
                complex(el[2]).imag,
                el[-1],
            ]
            ret2[indice] = valeur.copy()
        #print(ret2[0][1]==0 and ret2[0][3]==0)
        return ret2
    else:
        indices = [index for index, value in enumerate(L) if value == True]

        bord[1][1] += sqrt(2) / 4 * (-1) ** indices[0] * 1 / dx
        bord[1][2] += (sqrt(2) / 4 * (-1) ** indices[0] * 1 / dx)*complex(0,1)

        bord[2][1] += -sqrt(2) / 4 * (-1) ** indices[0] * 1 / dx
        bord[2][2] += (-sqrt(2) / 4 * (-1) ** indices[0] * 1 / dx)*complex(0,1)

        bord[3][1] += sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy
        bord[3][2] += (sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy)*complex(0,1)

        bord[4][1] += -sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy
        bord[4][2] += (-sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy)*complex(0,1)

        # réinjection et suppression du point fantôme

        for i in range(5):
            if i - 1 not in indices:
                ret[i][1] += -bord[i][1] * (
                    ret[indices[0] + 1][1] / bord[indices[0] + 1][1]
                    + ret[indices[0] + 1][2] / bord[indices[0] + 1][2]
                )
                ret[i][2] += -bord[i][2] * (
                    ret[indices[0] + 1][1] / bord[indices[0] + 1][1]
                    + ret[indices[0] + 1][2] / bord[indices[0] + 1][2]
                )

                ret[i][1] += -bord[i][1] * (
                    ret[indices[1] + 1][1] / bord[indices[1] + 1][1]
                    + ret[indices[1] + 1][2] / bord[indices[1] + 1][2]
                )
                ret[i][2] += -bord[i][2] * (
                    ret[indices[1] + 1][1] / bord[indices[1] + 1][1]
                    + ret[indices[1] + 1][2] / bord[indices[1] + 1][2]
                )
        del ret[indices[0] + 1]  # supprime le point fantôme de la liste ret
        del ret[indices[1]]  # supprime le point fantôme de la liste ret

        for indice in range(len(ret)):
            el = ret[indice]
            valeur = [
                el[0],
                complex(el[1]).real,
                complex(el[1]).imag,
                complex(el[2]).real,
                complex(el[2]).imag,
                el[-1],
            ]
            ret[indice] = valeur.copy()
        return ret


def BC_out(
    T, i, j, z=1
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_12
    dx = f_dx(z)
    dy = f_dy(z)
    L = [
        T[i, j + 1]["flag"] == "fantome",  # droite
        T[i, j - 1]["flag"] == "fantome",  # gauche
        T[i + 1, j]["flag"] == "fantome",  # bas
        T[i - 1, j]["flag"] == "fantome",  # haut
    ]
    l = sum(L)

    bord = [
        [
            (i, j),
            complex(0, 1) * k + 2 * complex(0, 1) * k0 * M0,
            (-complex(0, 1) * k + 2 * complex(0, 1) * k0 * M0) * complex(0, 1),
        ],
        [(i, j + 1), -(M0**2) / dx, -(M0**2) / dx],
        [(i, j - 1), M0**2 / dx, M0**2 / dx],
        [(i + 1, j), 0, 0],
        [(i - 1, j), 0, 0],
    ]  # donne l'équation des conditions sur le bord out

    D = ["in", "CGB", "CGH"]

    ret = [
        [
            (i, j),
            k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2,
            (k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2) * complex(0, 1),
            T[i, j]["flag"] == "frontiere" and T[i, j] in D,
        ],
        [
            (i, j + 1),
            1 / dx**2 - (M0 / dx) ** 2 + complex(0, 1) * k0 * (M0 / dx),
            (1 / dx**2 - (M0 / dx) ** 2 + complex(0, 1) * k0 * (M0 / dx))
            * complex(0, 1),
            T[i, j]["flag"] == "frontiere" and T[i, j + 1] in D,
        ],
        [
            (i, j - 1),
            1 / dx**2 - (M0 / dx) ** 2 - complex(0, 1) * k0 * (M0 / dx),
            (1 / dx**2 - (M0 / dx) ** 2 - complex(0, 1) * k0 * (M0 / dx))
            * complex(0, 1),
            T[i, j]["flag"] == "frontiere" and T[i, j - 1] in D,
        ],
        [(i + 1, j), 0, 0, T[i, j]["flag"] == "frontiere" and T[i + 1, j] in D],
        [(i - 1, j), 0, 0, T[i, j]["flag"] == "frontiere" and T[i - 1, j] in D],
    ] 

    if l == 1:
        r = L.index(True)

        # ajout de la dérivée sur n

        if r < 2:
            bord[1][1] += 1 / 2 * (-1) ** r * 1 / dx
            bord[1][2] += (1 / 2 * (-1) ** r * 1 / dx)*complex(0,1)

            bord[2][1] += -1 / 2 * (-1) ** r * 1 / dx
            bord[2][2] += (-1 / 2 * (-1) ** r * 1 / dx)*complex(0,1)
        else:
            bord[3][1] += -1 / 2 * (-1) ** r * 1 / dy
            bord[3][2] += (-1 / 2 * (-1) ** r * 1 / dy)*complex(0,1)

            bord[4][1] += 1 / 2 * (-1) ** r * 1 / dy
            bord[4][2] += (1 / 2 * (-1) ** r * 1 / dy)*complex(0,1)

        # réinjection et suppression du point fantôme

        for i in range(5):
            if i != r + 1:
                A=np.array([[bord[r+1][1].real,bord[r+1][2].real],[bord[r+1][1].imag,bord[r+1][2].imag]])
                B=-np.array([[bord[i][1].real,bord[i][2].real],[bord[i][1].imag,bord[i][2].imag]])
                C=np.linalg.inv(A)@B
                ret[i][1]+=ret[r+1][1]*C[0,0]+ret[r+1][2]*C[1,0]
                ret[i][2]+=ret[r+1][1]*C[0,1]+ret[r+1][2]*C[1,1]
                # ret[i][1] += -bord[i][1] * (
                #     ret[r + 1][1] / bord[r + 1][1] + ret[r + 1][2] / bord[r + 1][2]
                # )  # à vérifier
                # ret[i][2] += -bord[i][2] * (
                #     ret[r + 1][1] / bord[r + 1][1] + ret[r + 1][2] / bord[r + 1][2]
                # )  # à vérifier

        del ret[r + 1]  # supprime le point fantôme de la liste ret

        ret2 = ret.copy()
        for indice in range(len(ret)):
            el = ret[indice]
            valeur = [
                el[0],
                complex(el[1]).real,
                complex(el[1]).imag,
                complex(el[2]).real,
                complex(el[2]).imag,
                el[-1],
            ]
            ret2[indice] = valeur.copy()
        #print(ret2[0][1]==0 and ret2[0][3]==0)
        return ret2
    else:
        indices = [index for index, value in enumerate(L) if value == True]

        bord[1][1] += sqrt(2) / 4 * (-1) ** indices[0] * 1 / dx
        bord[1][2] += (sqrt(2) / 4 * (-1) ** indices[0] * 1 / dx)*complex(0,1)

        bord[2][1] += -sqrt(2) / 4 * (-1) ** indices[0] * 1 / dx
        bord[2][2] += (-sqrt(2) / 4 * (-1) ** indices[0] * 1 / dx)*complex(0,1)

        bord[3][1] += -sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy
        bord[3][2] += (-sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy)*complex(0,1)

        bord[4][1] += sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy
        bord[4][2] += (sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy)*complex(0,1)

        # réinjection et suppression du point fantôme

        for i in range(5):
            if i - 1 not in indices:
                ret[i][1] += -bord[i][1] * (
                    ret[indices[0] + 1][1] / bord[indices[0] + 1][1]
                    + ret[indices[0] + 1][2] / bord[indices[0] + 1][2]
                )  # à vérifier
                ret[i][2] += -bord[i][2] * (
                    ret[indices[0] + 1][1] / bord[indices[0] + 1][1]
                    + ret[indices[0] + 1][2] / bord[indices[0] + 1][2]
                )  # à vérifier

                ret[i][1] += -bord[i][1] * (
                    ret[indices[1] + 1][1] / bord[indices[1] + 1][1]
                    + ret[indices[1] + 1][2] / bord[indices[1] + 1][2]
                )  # à vérifier
                ret[i][2] += -bord[i][2] * (
                    ret[indices[1] + 1][1] / bord[indices[1] + 1][1]
                    + ret[indices[1] + 1][2] / bord[indices[1] + 1][2]
                )  # à vérifier
        del ret[indices[0] + 1]  # supprime le point fantôme de la liste ret
        del ret[indices[1]]  # supprime le point fantôme de la liste ret

        for indice in range(len(ret)):
            el = ret[indice]
            valeur = [
                el[0],
                complex(el[1]).real,
                complex(el[1]).imag,
                complex(el[2]).real,
                complex(el[2]).imag,
                el[-1],
            ]
            ret[indice] = valeur.copy()
        #print(ret[0][1]==0 and ret[0][3]==0)
        return ret


def omega(T, i, j, z=1):
    dx = f_dx(z)
    dy = f_dy(z)
    D = ["in", "CGB", "CGH"]

    L= [
        [
            (i, j),
            complex((k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2)).real,
            complex((k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2)).imag,
            complex(
                (k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2)
                * complex(0, 1)
            ).real,
            complex(
                (k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2)
                * complex(0, 1)
            ).imag,
            T[i, j]["flag"] == "frontiere" and T[i, j] in D,
        ],
        [
            (i, j + 1),
            complex(1 / dx**2 - (M0 / dx) ** 2 + complex(0, 1) * k0 * (M0 / dx)).real,
            complex(1 / dx**2 - (M0 / dx) ** 2 + complex(0, 1) * k0 * (M0 / dx)).imag,
            complex(
                (1 / dx**2 - (M0 / dx) ** 2 + complex(0, 1) * k0 * (M0 / dx))
                * complex(0, 1)
            ).real,
            complex(
                (1 / dx**2 - (M0 / dx) ** 2 + complex(0, 1) * k0 * (M0 / dx))
                * complex(0, 1)
            ).imag,
            T[i, j]["flag"] == "frontiere" and T[i, j + 1] in D,
        ],
        [
            (i, j - 1),
            complex(1 / dx**2 - (M0 / dx) ** 2 - complex(0, 1) * k0 * (M0 / dx)).real,
            complex(1 / dx**2 - (M0 / dx) ** 2 - complex(0, 1) * k0 * (M0 / dx)).imag,
            complex(
                (1 / dx**2 - (M0 / dx) ** 2 - complex(0, 1) * k0 * (M0 / dx))
                * complex(0, 1)
            ).real,
            complex(
                (1 / dx**2 - (M0 / dx) ** 2 - complex(0, 1) * k0 * (M0 / dx))
                * complex(0, 1)
            ).imag,
            T[i, j]["flag"] == "frontiere" and T[i, j - 1] in D,
        ],
        [(i + 1, j), 0, 0, 0, 0, T[i, j]["flag"] == "frontiere" and T[i + 1, j] in D],
        [(i - 1, j), 0, 0, 0, 0, T[i, j]["flag"] == "frontiere" and T[i - 1, j] in D],
    ]
    return L


def BC_in(
    T, i, j
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_in
    return [[(i, j), 1, 0, 0, 1, False]]


def BC_coin(
    T, i, j, chi, z=1
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_12
    dx = f_dx(z)
    dy = f_dy(z)
    L = [
        T[i, j + 1]["flag"] == "fantome",  # droite
        T[i, j - 1]["flag"] == "fantome",  # gauche
        T[i + 1, j]["flag"] == "fantome",  # bas
        T[i - 1, j]["flag"] == "fantome",  # haut
    ]

    indices = [index for index, value in enumerate(L) if value == True]

    if 1 in indices:  # cas où on a un coin à sur le bord in à gauche
        return [[(i, j), 1, 0, 0, 1, False]]

    # sinon on est dans le cas d'un coin entre 1,2 et out donc haut ou bas à droite

    bord_12 = [
        [
            (i, j),
            -complex(0, 1) * chi * alphaR
            + chi * alphaI
            - 2 * complex(0, 1) * chi * alphaR * (M0 / (k0 * dx)) ** 2
            + 2 * chi * alphaI * (M0 / (k0 * dx)) ** 2
            - 2 * complex(0, 1) * k0 * M0,
            (-complex(0, 1) * chi * alphaR
            + chi * alphaI
            + 2 * complex(0, 1) * chi * alphaR * (M0 / (k0 * dx)) ** 2
            - 2 * chi * alphaI * (M0 / (k0 * dx)) ** 2
            - 2 * complex(0, 1) * k0 * M0) * complex(0, 1),
        ],
        [
            (i, j + 1),
            (2 * chi * alphaR * (M0 / k0) + M0**2+2*chi*alphaI*M0/k0) / (2 * dx)
            + (complex(0, 1) * alphaR - alphaI) * chi * M0**2 / (k0 * dx) ** 2
            + M0**2 / (2 * dx),
            ((2 * chi * alphaR * (M0 / k0) + M0**2+2*chi*alphaI*M0/k0) / (2 * dx)
            - (complex(0, 1) * alphaR - alphaI) * chi * M0**2 / (k0 * dx) ** 2
            + M0**2 / (2 * dx)) * complex(0, 1),
        ],
        [
            (i, j - 1),
            -(
                (2 * chi * alphaR * (M0 / k0) + M0**2+2*chi*alphaI*M0/k0) / (2 * dx)
                + (complex(0, 1) * alphaR - alphaI) * chi * M0**2 / (k0 * dx) ** 2
                + M0**2 / (2 * dx)
            ),
            -(
                (2 * chi * alphaR * (M0 / k0) + M0**2+2*chi*alphaI*M0/k0) / (2 * dx)
                - (complex(0, 1) * alphaR - alphaI) * chi * M0**2 / (k0 * dx) ** 2
                + M0**2 / (2 * dx)
            ) * complex(0, 1),
        ],
        [(i + 1, j), 0, 0],
        [(i - 1, j), 0, 0],
    ]  # donne l'équation des conditions sur le bord 1 ou 2

    bord_out = [
        [
            (i, j),
            complex(0, 1) * k + 2 * complex(0, 1) * k0 * M0,
            (-complex(0, 1) * k + 2 * complex(0, 1) * k0 * M0) * complex(0, 1),
        ],
        [(i, j + 1), -(M0**2) / dx, (-(M0**2) / dx) * complex(0, 1)],
        [(i, j - 1), M0**2 / dx, (M0**2 / dx) * complex(0, 1)],
        [(i + 1, j), 0, 0],
        [(i - 1, j), 0, 0],
    ]  # donne l'équation des conditions sur le bord out

    D = ["in", "CGB", "CGH"]

    ret = [
        [
            (i, j),
            k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2,
            (k0**2 - 2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2) * complex(0, 1),
            T[i, j]["flag"] == "frontiere" and T[i, j] in D,
        ],
        [
            (i, j + 1),
            1 / dx**2 - (M0 / dx) ** 2 + complex(0, 1) * k0 * (M0 / dx),
            (1 / dx**2 - (M0 / dx) ** 2 + complex(0, 1) * k0 * (M0 / dx))
            * complex(0, 1),
            T[i, j]["flag"] == "frontiere" and T[i, j + 1] in D,
        ],
        [
            (i, j - 1),
            1 / dx**2 - (M0 / dx) ** 2 - complex(0, 1) * k0 * (M0 / dx),
            (1 / dx**2 - (M0 / dx) ** 2 - complex(0, 1) * k0 * (M0 / dx))
            * complex(0, 1),
            T[i, j]["flag"] == "frontiere" and T[i, j - 1] in D,
        ],
        [(i + 1, j), 0, 0, T[i, j]["flag"] == "frontiere" and T[i + 1, j] in D],
        [(i - 1, j), 0, 0, T[i, j]["flag"] == "frontiere" and T[i - 1, j] in D],
    ] 

    # ajout de la dérivée sur n

    bord_out[1][1] += sqrt(2) / 4 * 1 / dx
    bord_out[1][2] += (sqrt(2) / 4 * 1 / dx)*complex(0,1)

    bord_out[2][1] += -sqrt(2) / 4 * 1 / dx
    bord_out[2][2] += (-sqrt(2) / 4 * 1 / dx)*complex(0,1)

    bord_12[3][1] += sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy
    bord_12[3][2] += (sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy)*complex(0,1)

    bord_12[4][1] += -sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy
    bord_12[4][2] += (-sqrt(2) / 4 * (-1) ** indices[1] * 1 / dy)*complex(0,1)

    # réinjection et suppression du point fantôme

    for i in range(5):
        if i - 1 not in indices:
            ret[i][1] += -bord_out[i][1] * (
                ret[indices[0] + 1][1] / bord_out[indices[0] + 1][1]
                + ret[indices[0] + 1][2] / bord_out[indices[0] + 1][2]
            )  # à vérifier
            ret[i][2] += -bord_out[i][2] * (
                ret[indices[0] + 1][1] / bord_out[indices[0] + 1][1]
                + ret[indices[0] + 1][2] / bord_out[indices[0] + 1][2]
            )  # à vérifier

            ret[i][1] += -bord_12[i][1] * (
                ret[indices[1] + 1][1] / bord_12[indices[1] + 1][1]
                + ret[indices[1] + 1][2] / bord_12[indices[1] + 1][2]
            )  # à vérifier
            ret[i][2] += -bord_12[i][2] * (
                ret[indices[1] + 1][1] / bord_12[indices[1] + 1][1]
                + ret[indices[1] + 1][2] / bord_12[indices[1] + 1][2]
            )  # à vérifier

    del ret[indices[0] + 1]  # supprime le point fantôme de la liste ret
    del ret[indices[1]]  # supprime le point fantôme de la liste ret

    for indice in range(len(ret)):
        el = ret[indice]
        valeur = [
            el[0],
            complex(el[1]).real,
            complex(el[1]).imag,
            complex(el[2]).real,
            complex(el[2]).imag,
            el[-1],
        ]
        ret[indice] = valeur.copy()
    #print(ret[0][1]==0 and ret[0][3]==0)
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
                index for index, sublist in enumerate(L) if sublist[5] == True
            ]  # donne les indices des voisins de i,j qui sont dans Gamma_in
            a = 0  # a est un compteur pour savoir si B[indice équivalent à la ligne i,j] existe déjà
            for r in range(
                0, len(L)
            ):  # L[r][0] parcourt les indices des voisins de i,j à prendre en compte
                (u, v) = L[r][0]
                q = vectorise_adj(u, v, p)
                if r in K:  # on voit si le voisin (u,v) est dans Gamma_in
                    if a == 0:  # B[equivalent i,j] n'existe pas encore
                        B.append(b[k] - L[r][1] * b[q] - L[r][3] * b[q])# on fait le relevé de dirichlet
                    else:
                        B[-1] += -L[r][1] * b[q] - L[r][3] * b[q]
                        # B[...] existe déjà et on lui ajoute la partie correspondante à un autre voisin (u,v)
                    a += 1
                else:  # si (u,v) n'est pas dans Gamma_in on ajoute sa valeur correspondante dans la matrice
                    A.append(L[r][1])
                    A.append(L[r][3])
                    JA.append(q)
                    JA.append(q + 1)
            IA.append(len(A))
            
            for r in range(
                0, len(L)
            ):  # L[r][0] parcourt les indices des voisins de i,j à prendre en compte
                (u, v) = L[r][0]
                q = vectorise_adj(u, v, p)
                if r in K:  # on voit si le voisin (u,v) est dans Gamma_in
                    if a == 0:  # B[equivalent i,j] n'existe pas encore
                        B.append(b[k + 1] - L[r][2] * b[q + 1] - L[r][4] * b[q + 1])  # on fait le relevé de dirichlet
                    else:
                        B[-1] += -L[r][2] * b[q + 1] - L[r][4] * b[q + 1]
                        # B[...] existe déjà et on lui ajoute la partie correspondante à un autre voisin (u,v)
                    a += 1
                else:  # si (u,v) n'est pas dans Gamma_in on ajoute sa valeur correspondante dans la matrice
                    A.append(L[r][2])
                    A.append(L[r][4])
                    JA.append(q)
                    JA.append(q + 1)
            IA.append(len(A))
            
            if (a == 0):  # si aucun voisin n'est dans Gamma_in B[...] on lui donne la valeur de b[k]
                B.append(b[k])
                B.append(b[k + 1])
            IB.append((i, j))

        elif T[i, j]["flag"] == "frontiere":
            if T[i, j]["BC"] in ["in", "CGH", "CGB"]:
                L = BC_in(T, i, j)
            if T[i, j]["BC"] == "out":
                L = BC_out(T, i, j, z)
            if T[i, j]["BC"] == "12":
                L = BC_12(T, i, j, chi(i, j), z)
            if T[i, j]["BC"] == "Neumann":
                L = BC_12(T, i, j, chi(i, j), z)
            if T[i, j]["BC"] in ["CDB", "CDH"]:
                L = BC_coin(T, i, j, chi(i, j), z)
            K = [index for index, sublist in enumerate(L) if sublist[2] == True]
            a = 0
            for r in range(
                0, len(L)
            ):  # L[r][0] parcourt les indices des voisins de i,j à prendre en compte
                (u, v) = L[r][0]
                q = vectorise_adj(u, v, p)
                if r in K:  # on voit si le voisin (u,v) est dans Gamma_in
                    if a == 0:  # B[equivalent i,j] n'existe pas encore
                        B.append(b[k] - L[r][1] * b[q] - L[r][3] * b[q])# on fait le relevé de dirichlet
                    else:
                        B[-1] += -L[r][1] * b[q] - L[r][3] * b[q]
                        # B[...] existe déjà et on lui ajoute la partie correspondante à un autre voisin (u,v)
                    a += 1
                else:  # si (u,v) n'est pas dans Gamma_in on ajoute sa valeur correspondante dans la matrice
                    A.append(L[r][1])
                    A.append(L[r][3])
                    JA.append(q)
                    JA.append(q + 1)
            IA.append(len(A))
            
            for r in range(
                0, len(L)
            ):  # L[r][0] parcourt les indices des voisins de i,j à prendre en compte
                (u, v) = L[r][0]
                q = vectorise_adj(u, v, p)
                if r in K:  # on voit si le voisin (u,v) est dans Gamma_in
                    if a == 0:  # B[equivalent i,j] n'existe pas encore
                        B.append(b[k + 1] - L[r][2] * b[q + 1] - L[r][4] * b[q + 1])  # on fait le relevé de dirichlet
                    else:
                        B[-1] += -L[r][2] * b[q + 1] - L[r][4] * b[q + 1]
                        # B[...] existe déjà et on lui ajoute la partie correspondante à un autre voisin (u,v)
                    a += 1
                else:  # si (u,v) n'est pas dans Gamma_in on ajoute sa valeur correspondante dans la matrice
                    A.append(L[r][2])
                    A.append(L[r][4])
                    JA.append(q)
                    JA.append(q + 1)
            IA.append(len(A))
            
            if (a == 0):  # si aucun voisin n'est dans Gamma_in B[...] on lui donne la valeur de b[k]
                B.append(b[k])
                B.append(b[k + 1])
            IB.append((i, j))
        else:
            IA.append(len(A))
            IA.append(len(A))
        #print(k)
        k += 2  # A chaque itération on coupe la partie imaginaire et la partie réelle en 2 valeurs distinctes
    # Conversion des listes en tableaux NumPy
    A = np.array(A, dtype=complex)
    IA = np.array(IA, dtype=int)
    JA = np.array(JA, dtype=int)

    # Crée une matrice CSR
    csrmatrix = scipy.sparse.csr_matrix((A, JA, IA), shape=(2 * n * p, 2 * n * p))

    # Supprimer les lignes et les colonnes nulles
    non_zero_rows = csrmatrix.getnnz(axis=1) > 0
    non_zero_cols = csrmatrix.getnnz(axis=0) > 0
    csrmatrix = csrmatrix[non_zero_rows][:, non_zero_cols]

    return (csrmatrix, B, IB)


def inverse_csr(A, b):
    return spsolve(A, b)  # b doit surement être sous format csr


matrice_reacteur = Mapping.mapping_BC(
    Mapping.color_to_flag(Mapping.mapping(Mapping.png_to_rgb_matrix("exemple7.png")))
)



                




def second_membre(u,T, g, z=1):  # check
    dx = f_dx(z)
    dy = f_dy(z)
    D = ["in", "CGH", "CGB"]
    #matrice_reacteur = Mapping.mapping_BC(Mapping.color_to_flag(Mapping.mapping(Mapping.png_to_rgb_matrix(image))))
    matrice_reacteur=T
    n, p = matrice_reacteur.shape[0], matrice_reacteur.shape[1]
    b = np.zeros(2 * n * p, dtype=complex)
    k=0
    while k<2*n*p:
        (i,j)=matricise_adj(k,p)
        if (
            matrice_reacteur[i, j] is not None
            and "BC" in matrice_reacteur[i, j]
            and matrice_reacteur[i, j]["BC"] in D
        ):
            # à changer selon la discretisation de g
            b[vectorise_adj(i, j, p)] = 0
            b[vectorise_adj(i, j, p) + 1] = 0
        elif matrice_reacteur[i, j] is not None and matrice_reacteur[i, j]['flag'] in ['EDP','frontiere']:
            b[vectorise_adj(i, j, p)] = g(u,dx * j, dy * (n - i),T, z)[0]
            b[vectorise_adj(i, j, p) + 1] = g(u,dx * j, dy * (n - i),T, z)[1]
        k+=2
    return b


def solution(u,
    T, g, chi, z=1
):  # image est l'image du reacteur et g la fonction valeur en gamma in
    # matrice_reacteur = Mapping.mapping_BC(
    # Mapping.color_to_flag(Mapping.mapping(Mapping.png_to_rgb_matrix(image))))
    matrice_reacteur=T
    n, p = matrice_reacteur.shape[0], matrice_reacteur.shape[1]
    sm=second_membre(u,matrice_reacteur, g, z)
    M=matrice_to_csr(matrice_reacteur,chi,sm)
    # u = inverse_csr(matrice_to_csr(Mapping.mapping_BC(Mapping.color_to_flag(Mapping.mapping(matrice_reacteur))),chi,second_membre(image, g,z),)[0],matrice_to_csr(Mapping.mapping_BC(Mapping.color_to_flag(Mapping.mapping(matrice_reacteur))),chi,second_membre(image, g,z),)[1],)
    # IB = matrice_to_csr(Mapping.mapping_BC(Mapping.color_to_flag(Mapping.mapping(matrice_reacteur))),chi,second_membre(image, g,z),)[2]
    u = inverse_csr(M[0],M[1])
    IB = M[2]
    matrice_pression = np.zeros((n, p), dtype=complex)
    k = 0
    while k < len(
        IB
    ):  # sous hypothèse que ca marche bien et que tout élement non nul a sa partie imaginaire et sa partie réelle non nulles, sinon à modifier peut etre, en rajoutant dans IB si c'est une partie réelle ou une partie imaginaire
        (i, j) = IB[k]
        matrice_pression[i, j] = u[k] + complex(0, 1) * u[k + 1]
        k += 2
    #background = plt.imread(image)
    # # Créez un tracé de la matrice de pression avec une échelle de couleur
    # min_real = np.min(matrice_pression.real)
    # max_real = np.max(matrice_pression.real)
    # min_imag = np.min(matrice_pression.imag)
    # max_imag = np.max(matrice_pression.imag)
    # fig, (ax1, ax2) = plt.subplots(1, 2)

    # # Afficher la partie réelle dans le premier sous-graphique
    # im1 = ax1.imshow(matrice_pression.real, cmap="coolwarm", interpolation="none")
    # ax1.set_title("Partie Réelle")
    # im1.set_clim(vmin=min_real, vmax=max_real)  # Ajuster l'échelle de couleur
    # fig.colorbar(im1, ax=ax1)

    # # Afficher la partie imaginaire dans le deuxième sous-graphique
    # im2 = ax2.imshow(matrice_pression.imag, cmap="coolwarm", interpolation="none")
    # ax2.set_title("Partie Imaginaire")
    # im2.set_clim(vmin=min_imag, vmax=max_imag)  # Ajuster l'échelle de couleur
    # fig.colorbar(im2, ax=ax2)

    # # Ajuster l'espacement entre les sous-graphiques
    # plt.tight_layout()

    return matrice_pression


def d_adjoint_dx(u,T, g, chi, z=1):
    dx = f_dx(z)
    matrice_reacteur=T
    #matrice_reacteur = Mapping.png_to_rgb_matrix(image)
    n, p = matrice_reacteur.shape[0], matrice_reacteur.shape[1]
    u = solution(u,T, g, chi, z)
    D = np.zeros((n, p), dtype=complex)
    for k in range(0, n * p):
        (i, j) = Mapping.matricise(k, p)
        if j == 0:
            D[i, j] = (u[i, j + 1] - u[i, j]) / dx
        elif j == p - 1:
            D[i, j] = (u[i, j] - u[i, j - 1]) / dx
        else:
            D[i, j] = (u[i, j + 1] - u[i, j - 1]) / (2 * dx)
    return D


def chi(i, j):
    return 0


def fct_test(x, y,T, z=1):
#return -2p_conj
    # matrice_reacteur = Mapping.mapping_BC(
    # Mapping.color_to_flag(Mapping.mapping(Mapping.png_to_rgb_matrix(image))))
    matrice_reacteur=T
    x0, x1, y0, y1 = +inf, -inf, +inf, -inf

    for i in range(matrice_reacteur.shape[0]):
        for j in range(matrice_reacteur.shape[1]):
            if "BC" in matrice_reacteur[i,j] and matrice_reacteur[i,j]["BC"] == "in":
                if j < x0:
                    x0 = j
                if i < y0:
                    y0 = i
                if i > y1:
                    y1 = i
            elif "BC" in matrice_reacteur[i,j] and matrice_reacteur[i,j]["BC"] == "out":
                if j > x1:
                    x1 = j
                if i < y0:
                    y0 = i
                if i > y1:
                    y1 = i
              
    x0 = x0 * f_dx(z)
    x1 = x1 * f_dx(z)
    y0 = (y0-1) * f_dy(z)
    y1 = (y1+1) * f_dy(z)
    
    P = (
        -1/2 * (y-y0)**2 * (y-y1)**2 * (k0**2 * (x-x0)**2 * (x-x1) + 2*(1-M0**2) * ((x-x1) + 2*(x-x0)))
        - (x-x0)**2 * (x-x1) * ((y-y1)**2 + (y-y0)**2 + 4 * (y-y0)*(y-y1)),
        k0*M0*(y-y0)**2*(y-y1)**2 * (2*(x-x0)*(x-x1) + (x-x0)**2)
    )
    
    return (-2 * P[0], 2 * P[1])


# solution(matrice_reacteur, fct_test, chi, 1)
# plt.show()
# for i in range(8,23):
#     print('i=',i,'k=',vectorise_adj(i,8,120),u[i,8].real, second_membre("exemple7.png", fct_test,1)[vectorise_adj(i,8,120)],second_membre("exemple7.png", fct_test,1)[vectorise_adj(i,8,120)+1])