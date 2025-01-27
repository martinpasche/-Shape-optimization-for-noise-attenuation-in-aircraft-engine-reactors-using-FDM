from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import numpy as np
import Mapping
import importlib
from math import *
import scipy.sparse
import cmath
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import fonction_test
from parameters import *


def f_dx(k):
    return 0.01 / k


def f_dy(k):
    return 0.01 / k





# n = 10
# p = 10
def f_d(k):
    return 4 * f_dy(k)


# ATTENTION


def BC_out(
    T, i, j, z=1
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_out
    dx = f_dx(z)
    dy = f_dy(z)
    L = [
        T[i, j + 1]["flag"] == "fantome",
        T[i, j - 1]["flag"] == "fantome",
        T[i + 1, j]["flag"] == "fantome",
        T[i - 1, j]["flag"] == "fantome",
    ]
    D = ["in", "CGB", "CGH"]
    l = sum(L)

    if l == 1:
        r = L.index(True)

        if r < 2:
            return [
                [
                    (i, j),
                    complex(
                        -2 / dx**2
                        - 2 / dy**2
                        + k0**2
                        + 2 * (M0 / dx) ** 2
                        - 2 * M0 * k * k0 * (-1) ** r,
                        -(2 * k / dx) + ((M0**2 * 2 * k) / dx),
                    ),
                    T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] == "in",
                ],
                [
                    (i, j - (-1) ** r),
                    2 / dx**2 - 2 * M0**2 / dx**2,
                    T[i, j - (-1) ** r]["flag"] == "frontiere"
                    and T[i, j - (-1) ** r]["BC"] == "in",
                ],
                [
                    (i + 1, j),
                    1 / dy**2,
                    T[i + 1, j]["flag"] == "frontiere" and T[i + 1, j]["BC"] in D,
                ],
                [
                    (i - 1, j),
                    1 / dy**2,
                    T[i - 1, j]["flag"] == "frontiere" and T[i - 1, j]["BC"] in D,
                ],
            ]

        else:
            return [
                [
                    (i, j),
                    complex(
                        -2 / dx**2 - 2 / dy**2 + k0**2 + 2 * (M0 / dx) ** 2,
                        -(2 * k / dy),
                    ),
                    T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] == "in",
                ],
                [
                    (i, j + 1),
                    complex(1 / dx**2 - M0**2 / dx**2, -(k0 * M0 / dx)),
                    T[i, j + 1]["flag"] == "frontiere" and T[i, j + 1]["BC"] in D,
                ],
                [
                    (i, j - 1),
                    complex(1 / dx**2 - M0**2 / dx**2, +(k0 * M0 / dx)),
                    T[i, j - 1]["flag"] == "frontiere" and T[i, j - 1]["BC"] in D,
                ],
                [
                    (i - (-1) ** r, j),
                    2 / dy**2,
                    T[i - (-1) ** r, j]["flag"] == "frontiere"
                    and T[i - (-1) ** r, j]["BC"] in D,
                ],
            ]

    else:
        indices = [index for index, value in enumerate(L) if value == True]
        r = indices[0]
        s = indices[1]
        return [
            [
                (i, j),
                -2 / dx**2
                - 2 / dy**2
                + k0**2
                + 2 * (M0 / dx) ** 2
                - complex(0, 2 * dx * k)
                * complex((1 - M0**2) / dx**2, -((-1) ** r) * M0 * k0 / dx)
                - complex(0, 2 * k / dy),
                T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] in D,
            ],
            [
                (i - (-1) ** s, j),
                2 / dy**2,
                T[i - (-1) ** r, j]["flag"] == "frontiere"
                and T[i - (-1) ** r, j]["BC"] in D,
            ],
            [
                (i, j - (-1) ** r),
                2 / dx**2 - 2 * M0**2 / dx**2,
                T[i, j - (-1) ** r]["flag"] == "frontiere"
                and T[i, j - (-1) ** r]["BC"] in D,
            ],
        ]


def BC_12(
    T, i, j, chi, z=1
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_12
    dx = f_dx(z)
    dy = f_dy(z)
    L = [
        T[i, j + 1]["flag"] == "fantome",
        T[i, j - 1]["flag"] == "fantome",
        T[i + 1, j]["flag"] == "fantome",
        T[i - 1, j]["flag"] == "fantome",
    ]
    D = ["in", "CGB", "CGH"]
    l = sum(L)

    if l == 1:
        r = L.index(True)

        if r < 2:
            return [
                [
                    (i, j),
                    -2 / dx**2
                    - 2 / dy**2
                    + k0**2
                    + 2 * (M0 / dx) ** 2
                    - complex(1 / dx**2 - (M0 / dx) ** 2, -((-1) ** r) * M0 * k0 / dx)
                    * complex(0, chi * Z0 * (k0 + 2 * M0**2 / (k0 * dx**2)) / Z)
                    / complex(
                        0.5 / dx + ((-1) ** r) * Z0 * M0 * chi / (Z * dx),
                        -chi * Z0 * M0**2 / (Z * k0 * dx**2),
                    ),
                    T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] in D,
                ],
                [
                    (i, j - (-1) ** r),
                    complex(1 / dx**2 - M0**2 / dx**2, ((-1) ** r) * M0 * k0 / dx)
                    + complex(1 / dx**2 - (M0 / dx) ** 2, -((-1) ** r) * M0 * k0 / dx)
                    * complex(
                        0.5 / dx + ((-1) ** r) * Z0 * M0 * chi / (Z * dx),
                        chi * Z0 * M0**2 / (Z * k0 * dx**2),
                    )
                    / complex(
                        0.5 / dx + ((-1) ** r) * Z0 * M0 * chi / (Z * dx),
                        -chi * Z0 * M0**2 / (Z * k0 * dx**2),
                    ),
                    T[i, j - (-1) ** r]["flag"] == "frontiere"
                    and T[i, j - (-1) ** r]["BC"] in D,
                ],
                [
                    (i + 1, j),
                    1 / dy**2,
                    T[i + 1, j]["flag"] == "frontiere" and T[i + 1, j]["BC"] in D,
                ],
                [
                    (i - 1, j),
                    1 / dy**2,
                    T[i - 1, j]["flag"] == "frontiere" and T[i - 1, j]["BC"] in D,
                ],
            ]
        else:
            return [
                [
                    (i, j),
                    complex(
                        -2 / dx**2
                        - 2 / dy**2
                        + k0**2
                        + 2 * (M0 / dx) ** 2
                        - 4 * chi * Z0 * M0**2 / (Z * k0 * dy * dx**2),
                        -2 * chi * k0 * Z0 / (Z * dy),
                    ),
                    T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] in D,
                ],
                [
                    (i, j + 1),
                    complex(
                        1 / dx**2
                        - M0**2 / dx**2
                        - 2 * Z0 * M0 * chi / (dy * dx * Z),
                        -(k0 * M0 / dx)
                        + 2 * chi * Z0 * M0**2 / (dy * Z * k0 * dx**2),
                    ),
                    T[i, j + 1]["flag"] == "frontiere" and T[i, j + 1]["BC"] in D,
                ],
                [
                    (i, j - 1),
                    complex(
                        1 / dx**2
                        - M0**2 / dx**2
                        + 2 * Z0 * M0 * chi / (dy * dx * Z),
                        (k0 * M0 / dx)
                        + 2 * chi * Z0 * M0**2 / (dy * Z * k0 * dx**2),
                    ),
                    T[i, j - 1]["flag"] == "frontiere" and T[i, j - 1]["BC"] in D,
                ],
                [
                    (i - (-1) ** r, j),
                    2 / dy**2,
                    T[i - (-1) ** r, j]["flag"] == "frontiere"
                    and T[i - (-1) ** r, j]["BC"] in D,
                ],
            ]

    else:
        indices = [index for index, value in enumerate(L) if value == True]
        r = indices[0]
        s = indices[1]
        return [
            [
                (i, j),
                -2 / dx**2
                - 2 / dy**2
                + k0**2
                + 2 * (M0 / dx) ** 2
                + (
                    complex((1 - M0**2) / dx**2, -((-1) ** r) * M0 * k0 / dx)
                    + complex(0, 2 * k0 * chi * Z0 / (Z * dy))
                    * complex((M0 / (k0 * dx)) ** 2, (-1) ** r * M0 / (k0 * dx))
                )
                * complex(0, -(1 + 2 * (M0 / (k0 * dx)) ** 2) * k0 * Z0 * chi / Z)
                / complex(
                    1 / 2 * dx + (-1) ** r * chi * Z0 * M0 / (Z * dx),
                    -chi * Z0 * M0**2 / (Z * k0 * dx**2),
                )
                - (1 + 2 * (M0 / (k0 * dx)) ** 2)
                * complex(0, 2 * k0 * chi * Z0 / (Z * dy)),
                T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] in D,
            ],
            [
                (i - (-1) ** s, j),
                2 / dy**2,
                T[i - (-1) ** r, j]["flag"] == "frontiere"
                and T[i - (-1) ** r, j]["BC"] in D,
            ],
            [
                (i, j - (-1) ** r),
                complex(1 / dx**2 - 1 * M0**2 / dx**2, (-1) ** r * M0 * k0 / dx)
                + (
                    complex((1 - M0**2) / dx**2, -((-1) ** r) * M0 * k0 / dx)
                    + complex(0, 2 * k0 * chi * Z0 / (Z * dy))
                    * complex((M0 / (k0 * dx)) ** 2, (-1) ** r * M0 / (k0 * dx))
                )
                * (
                    1 / (2 * dx)
                    - complex(-((M0 / (k0 * dx)) ** 2), (-1) ** r * M0 / (k0 * dx))
                    * complex(k0 * Z0 * chi / Z)
                )
                / complex(
                    1 / (2 * dx) + (-1) ** r * chi * Z0 * M0 / (Z * dx),
                    -chi * Z0 * M0**2 / (Z * k0 * dx**2),
                )
                - complex(-((M0 / (k0 * dx)) ** 2), (-1) ** r * M0 / (k0 * dx))
                * complex(0, 2 * k0 * chi * Z0 / (Z * dy)),
                T[i, j - (-1) ** r]["flag"] == "frontiere"
                and T[i, j - (-1) ** r]["BC"] in D,
            ],
        ]


def BC_in(
    T, i, j, z=1
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur gamma_in
    return [[(i, j), 1, False]]


def BC_cdb(
    T, i, j, chi, z=1
):  # fonction pour le coin inférieur droit qui est à la fois sur gamma_out et gamma_12
    dx = f_dx(z)
    dy = f_dy(z)
    D = ["in", "CGB", "CGH"]
    if chi != 0:
        return [
            [
                (i, j),
                2 * (M0**2 - 1) / dx**2
                - 2 / dy**2
                + k0**2
                - complex(0, 2 * np.sqrt(2) * k / dy)
                + complex(
                    0,
                    chi * k0 * Z0 / Z + 2 * chi * Z0 * M0**2 / (Z * k0 * dx**2) - k,
                )
                * complex((1 - M0**2) / dx**2 - 1 / (dy * dx), -M0 * k0 / dx)
                / complex(
                    -chi * Z0 * M0 / (Z * dx), chi * Z0 * M0**2 / (Z * k0 * dx**2)
                ),
                T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] in D,
            ],
            [
                (i, j - 1),
                complex((1 - M0**2) / dx**2 + 1 / (dx * dy), M0 * k0 / dx)
                - complex(
                    chi * Z0 * M0 / (Z * dx), chi * Z0 * M0**2 / (Z * k0 * dx**2)
                )
                * complex((1 - M0**2) / dx**2 - 1 / (dx * dy), -M0 * k0 / dx)
                / complex(
                    -chi * Z0 * M0 / (Z * dx), chi * Z0 * M0**2 / (Z * k0 * dx**2)
                ),
                T[i, j - 1]["flag"] == "frontiere" and T[i, j - 1]["BC"] in D,
            ],
            [
                (i - 1, j),
                2 / dy**2,
                T[i + 1, j]["flag"] == "frontiere" and T[i + 1, j]["BC"] in D,
            ],
        ]
    else:
        return [
            [
                (i, j),
                2 * (M0**2 - 1) / dx**2
                - 2 / dy**2
                + k0**2
                - complex(0, 2 * np.sqrt(2) * k / dy),
                T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] in D,
            ],
            [
                (i, j - 1),
                complex((1 - M0**2) / dx**2 + 1 / (dx * dy), M0 * k0 / dx),
                T[i, j - 1]["flag"] == "frontiere" and T[i, j - 1]["BC"] in D,
            ],
            [
                (i - 1, j),
                2 / dy**2,
                T[i + 1, j]["flag"] == "frontiere" and T[i + 1, j]["BC"] in D,
            ],
        ]


def BC_cdh(
    T, i, j, chi, z=1
):  # fonction pour le coin supérieur droit qui est à la fois sur gamma_out et gamma_12
    dx = f_dx(z)
    dy = f_dy(z)
    D = ["in", "CGB", "CGH"]
    if chi != 0:
        return [
            [
                (i, j),
                2 * (M0**2 - 1) / dx**2
                - 2 / dy**2
                + k0**2
                - complex(0, 2 * np.sqrt(2) * k / dy)
                + complex(
                    0,
                    chi * k0 * Z0 / Z + 2 * chi * Z0 * M0**2 / (Z * k0 * dx**2) - k,
                )
                * complex((1 - M0**2) / dx**2 - 1 / (dy * dx), -M0 * k0 / dx)
                / complex(
                    -chi * Z0 * M0 / (Z * dx), chi * Z0 * M0**2 / (Z * k0 * dx**2)
                ),
                T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] in D,
            ],
            [
                (i, j - 1),
                complex((1 - M0**2) / dx**2 + 1 / (dx * dy), M0 * k0 / dx)
                - complex(
                    chi * Z0 * M0 / (Z * dx), chi * Z0 * M0**2 / (Z * k0 * dx**2)
                )
                * complex((1 - M0**2) / dx**2 - 1 / (dx * dy), -M0 * k0 / dx)
                / complex(
                    -chi * Z0 * M0 / (Z * dx), chi * Z0 * M0**2 / (Z * k0 * dx**2)
                ),
                T[i, j - 1]["flag"] == "frontiere" and T[i, j - 1]["BC"] in D,
            ],
            [
                (i + 1, j),
                2 / dy**2,
                T[i + 1, j]["flag"] == "frontiere" and T[i + 1, j]["BC"] in D,
            ],
        ]
    else:
        return [
            [
                (i, j),
                2 * (M0**2 - 1) / dx**2
                - 2 / dy**2
                + k0**2
                - complex(0, 2 * np.sqrt(2) * k / dy),
                T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] in D,
            ],
            [
                (i, j - 1),
                complex((1 - M0**2) / dx**2 + 1 / (dx * dy), M0 * k0 / dx),
                T[i, j - 1]["flag"] == "frontiere" and T[i, j - 1]["BC"] in D,
            ],
            [
                (i + 1, j),
                2 / dy**2,
                T[i + 1, j]["flag"] == "frontiere" and T[i + 1, j]["BC"] in D,
            ],
        ]


def BC_cgb(T, i, j, z=1):
    return [[(i, j), 1, False]]


def BC_cgh(T, i, j, z=1):
    return [[(i, j), 1, False]]


def omega(
    T, i, j, z=1
):  # donne les coordonnées et les valeurs associées à mettre dans la matrice lorsqu'on se situe sur omega
    dx = f_dx(z)
    dy = f_dy(z)
    D = ["in", "CGB", "CGH"]
    return [
        [
            (i, j),
            -2 / dx**2 - 2 / dy**2 + 2 * (M0 / dx) ** 2 + k0**2,
            T[i, j]["flag"] == "frontiere" and T[i, j]["BC"] in D,
        ],
        [
            (i, j + 1),
            complex(1 / dx**2 - (M0 / dx) ** 2, -M0 * k0 / dx),
            T[i, j + 1]["flag"] == "frontiere" and T[i, j + 1]["BC"] in D,
        ],
        [
            (i, j - 1),
            complex(1 / dx**2 - (M0 / dx) ** 2, M0 * k0 / dx),
            T[i, j - 1]["flag"] == "frontiere" and T[i, j - 1]["BC"] in D,
        ],
        [
            (i + 1, j),
            1 / dy**2,
            T[i + 1, j]["flag"] == "frontiere" and T[i + 1, j]["BC"] in D,
        ],
        [
            (i - 1, j),
            1 / dy**2,
            T[i - 1, j]["flag"] == "frontiere" and T[i - 1, j]["BC"] in D,
        ],
    ]


def chi(i, j):
    return 0


def G(i, j):
    return 1


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
        (i, j) = Mapping.matricise(k, p)
        if T[i, j] is not None and T[i, j]["flag"] == "EDP":
            L = omega(T, i, j, z)
            K = [
                index for index, sublist in enumerate(L) if sublist[2] == True
            ]  # donne les indices des voisins de i,j qui sont dans Gamma_in
            a = 0  # a est un compteur pour savoir si B[indice équivalent à la ligne i,j] existe déjà
            for r in range(
                0, len(L)
            ):  # L[r][0] parcourt les indices des voisins de i,j à prendre en compte
                (u, v) = L[r][0]
                q = Mapping.vectorise(u, v, p)
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

        if T[i, j] is not None and T[i, j]["flag"] == "frontiere":
            if T[i, j]["BC"] == "in":
                L = BC_in(T, i, j, z)
            if T[i, j]["BC"] == "out":
                L = BC_out(T, i, j, z)
            if T[i, j]["BC"] == "12":
                L = BC_12(T, i, j, chi(i, j), z)
            if T[i, j]["BC"] == "Neumann":
                L = BC_12(T, i, j, chi(i, j), z)
            if T[i, j]["BC"] == "CDB":
                L = BC_cdb(T, i, j, chi(i, j), z)
            if T[i, j]["BC"] == "CDH":
                L = BC_cdh(T, i, j, chi(i, j), z)
            if T[i, j]["BC"] == "CGH":
                L = BC_cgh(T, i, j, z)
            if T[i, j]["BC"] == "CGB":
                L = BC_cgb(T, i, j, z)
            K = [index for index, sublist in enumerate(L) if sublist[2] == True]
            a = 0
            for r in range(0, len(L)):
                (u, v) = L[r][0]
                q = Mapping.vectorise(u, v, p)
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

    return (csrmatrix, B, IB)


def inverse_csr(A, b):
    return spsolve(A, b)  # b doit surement être sous format csr


matrice_reacteur = Mapping.mapping_BC(
    Mapping.color_to_flag(Mapping.mapping(Mapping.png_to_rgb_matrix("exemple7.png")))
)


def second_membre(image, g, z=1):  # check
    dx = f_dx(z)
    dy = f_dy(z)
    D = ["in", "CGB", "CGH"]
    matrice_reacteur = Mapping.mapping_BC(Mapping.color_to_flag(Mapping.mapping(Mapping.png_to_rgb_matrix(image))))
    n, p = matrice_reacteur.shape[0], matrice_reacteur.shape[1]
    b = np.zeros(n * p, dtype=complex)
    for i in range(n):
        for j in range(p):
            if (
                matrice_reacteur[i, j] is not None
                and "BC" in matrice_reacteur[i, j]
                and matrice_reacteur[i, j]["BC"] in D
            ):
                # à changer selon la discretisation de g
                b[Mapping.vectorise(i, j, p)] = g(dx * j, dy * (n - i), "BC", z)
            else:
                b[Mapping.vectorise(i, j, p)] = g(dx * j, dy * (n - i), "gamma", z)
    return b


def fct_test(x, y, a, z=1):
    if a == "BC":
        return cos( pi * (y-8) / 14) * cmath.exp(complex(0, -k0 * x))
    else:
        return (
            cos(2 * pi * y / f_d(z))
            * cmath.exp(complex(0, -k * x))
            * (-k0 * k0 - 4 * pi * pi / f_d(z) ** 2 + k0 * k0 * (1 - M0) ** 2)
        )


def solution(
    image, g, chi, z=1
):  # image est l'image du reacteur et g la fonction valeur en gamma in
    matrice_reacteur = Mapping.png_to_rgb_matrix(image)
    n, p = matrice_reacteur.shape[0], matrice_reacteur.shape[1]
    u = inverse_csr(matrice_to_csr(Mapping.mapping_BC(Mapping.color_to_flag(Mapping.mapping(matrice_reacteur))),chi,second_membre(image, g,z),)[0],matrice_to_csr(Mapping.mapping_BC(Mapping.color_to_flag(Mapping.mapping(matrice_reacteur))),chi,second_membre(image, g,z),)[1],)
    IB = matrice_to_csr(Mapping.mapping_BC(Mapping.color_to_flag(Mapping.mapping(matrice_reacteur))),chi,second_membre(image, g,z),)[2]
    # u = inverse_csr(
    #     matrice_to_csr(
    #         matrice_reacteur,
    #         chi,
    #         second_membre(image, g, z),
    #     )[0],
    #     matrice_to_csr(
    #         matrice_reacteur,
    #         chi,
    #         second_membre(image, g, z),
    #     )[1],
    # )
    # IB = matrice_to_csr(
    #     matrice_reacteur,
    #     chi,
    #     second_membre(image, g, z),
    # )[2]
    matrice_pression = np.zeros((n, p), dtype=complex)
    for k in range(len(u)):
        (i, j) = IB[k]
        matrice_pression[i, j] = u[k]
    background = plt.imread(image)
    # Créez un tracé de la matrice de pression avec une échelle de couleur
    min_real = np.min(matrice_pression.real)
    max_real = np.max(matrice_pression.real)
    min_imag = np.min(matrice_pression.imag)
    max_imag = np.max(matrice_pression.imag)
    fig, (ax1, ax2) = plt.subplots(1, 2)

    # Afficher la partie réelle dans le premier sous-graphique
    im1 = ax1.imshow(matrice_pression.real, cmap="coolwarm", interpolation="none",origin="lower")
    ax1.set_title("Partie Réelle solution numérique")
    im1.set_clim(vmin=min_real, vmax=max_real)  # Ajuster l'échelle de couleur
    fig.colorbar(im1, ax=ax1)

    # Afficher la partie imaginaire dans le deuxième sous-graphique
    im2 = ax2.imshow(matrice_pression.imag, cmap="coolwarm", interpolation="none",origin="lower")
    ax2.set_title("Partie Imaginaire solution numérique")
    im2.set_clim(vmin=min_imag, vmax=max_imag)  # Ajuster l'échelle de couleur
    fig.colorbar(im2, ax=ax2)

    # # Ajuster l'espacement entre les sous-graphiques
    # plt.tight_layout()
    return matrice_pression, fig, ax1, ax2


def d_solution_dx(image, g, chi, z=1):
    dx = f_dx(z)
    # matrice_reacteur = Mapping.png_to_rgb_matrix(image)
    n, p = matrice_reacteur.shape[0], matrice_reacteur.shape[1]
    u = solution(image, g, chi, z)
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


def matrice_test_F(l):
    matrice = np.arange(1, l * l + 1).reshape((l, l))
    # Remplacer la ligne 5 par des zéros
    matrice[4, :] = 0

    # Remplacer la colonne 5 par des zéros
    matrice[:, 4] = 0

    print(matrice)

    # Supprimer les lignes et les colonnes nulles
    non_zero_rows = matrice.getnnz(axis=1) > 0
    non_zero_cols = matrice.getnnz(axis=0) > 0
    matrice = matrice[non_zero_rows][:, non_zero_cols]

    print(matrice)


def comparaison():
    return 0


# matrice_test_F(10)


# print(matrice_to_csr(matrice_reacteur,chi).get_shape())
# M = matrice_to_csr(matrice_reacteur, chi, second_membre("exemple7.png",fct_test))
# print(inverse_csr(M[0],M[1]))
#print(solution("exemple7.png",fct_test,chi))

mat, fig1, ax1, ax2 = solution("exemple7.png", fct_test, chi)
fig2, ax3, ax4 = fonction_test.verification()

# Ajouter des titres aux graphiques de la première figure
ax1.set_title("Partie Réelle calcué")
ax2.set_title("Partie Imaginaire calculé")

# Ajouter des titres aux graphiques de la deuxième figure
ax3.set_title("partie réelle analityque")
ax4.set_title("partie imaginaire analityque")


plt.show()
