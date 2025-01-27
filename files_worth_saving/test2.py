import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import code2


def norme_infini_ensemble(fonction, x_range, y_range, dx, dy):
    max_valeur = -float("inf")

    for x in np.arange(x_range[0], x_range[1], dx):
        for y in np.arange(y_range[0], y_range[1], dy):
            valeur = fonction(x, y)
            if abs(valeur) > max_valeur:
                max_valeur = abs(valeur)

    return max_valeur


def norme_L1_ensemble(fonction, x_range, y_range, dx, dy):
    somme = 0.0

    for x in np.arange(x_range[0], x_range[1], dx):
        for y in np.arange(y_range[0], y_range[1], dy):
            somme += abs(fonction(x, y))

    return somme


def norme_L2_ensemble(fonction, x_range, y_range, dx, dy):
    somme_carres = 0.0

    for x in np.arange(x_range[0], x_range[1], dx):
        for y in np.arange(y_range[0], y_range[1], dy):
            somme_carres += abs(fonction(x, y)) ** 2

    norme_L2 = np.sqrt(somme_carres)
    return norme_L2


def derivee_x(fonction, dx, x1, y1):
    return (fonction(x1 + dx, y1) - fonction(x1 - dx, y1)) / (2 * dx)


def derivee_y(fonction, dy, x1, y1):
    return (fonction(x1, y1 + dy) - fonction(x1, y1 - dy)) / (2 * dy)


def norme_H1_ensemble(fonction, x_range, y_range, dx, dy):
    somme_carres = 0.0
    somme_carres1 = 0.0
    for x in np.arange(x_range[0], x_range[1], dx):
        for y in np.arange(y_range[0], y_range[1], dy):
            somme_carres += (
                abs(derivee_x(fonction, dx, x, y)) ** 2
                + abs(derivee_y(fonction, dy, x, y)) ** 2
            )
            somme_carres1 += abs(fonction(x, y)) ** 2
    return np.sqrt(somme_carres1) + np.sqrt(somme_carres)


n = 120 * 2
p = 30 * 2
L = 105 # longueur exemple 7 en pixel
h = 15  # hauteur exemple 7 en pixel
k = code2.k


def verification(z=1):
    dx = code2.f_dx(z)
    dy = code2.f_dy(z)
    d = code2.f_d(z)
    # Créez des valeurs pour x et y
    x = np.linspace(0, n, n)
    y = np.linspace(0, p, p)
    X, Y = np.meshgrid(x, y)
    #print(X.shape)
    # Calculez votre fonction Z(x, y) - vous devez remplacer cette fonction par la vôtre
    # Exemple : une fonction simple
    Z = np.zeros((p, n), dtype=complex)
    condition = (X >= 8) & (X <= (112)) & (Y >= 8) & (Y <= 22)
    for i in range(0, p):
        for j in range(0, n):
            if condition[i, j] == True:
                Z[i, j] = np.cos(2 * np.pi / d * Y[i, j] * dy) * np.exp(
                    -complex(0, k * X[i, j] * dx)
                )
    # Obtenez la partie réelle et la partie imaginaire de Z
    Z_real = np.real(Z)
    Z_imag = np.imag(Z)
    # Créez une figure avec deux sous-graphiques côte à côte
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Graphique de la partie réelle
    im1 = ax1.imshow(
        Z_real,
        extent=(x.min(), x.max(), y.min(), y.max()),
        origin="lower",
        cmap="coolwarm",
    )
    ax1.set_title("Partie Réelle")
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    fig.colorbar(im1, ax=ax1)

    # Graphique de la partie imaginaire
    im2 = ax2.imshow(
        Z_imag,
        extent=(x.min(), x.max(), y.min(), y.max()),
        origin="lower",
        cmap="coolwarm",
    )
    ax2.set_title("Partie Imaginaire")
    ax2.set_xlabel("X")
    ax2.set_ylabel("Y")
    fig.colorbar(im2, ax=ax2)

    # Affichez la figure
    plt.show()


matrice_reacteur = code2.matrice_reacteur
chi = code2.chi


def P(x, y, z=1):
    dx = code2.f_dx(z)
    dy = code2.f_dy(z)
    M = code2.matrice_to_csr(
        matrice_reacteur, chi, code2.second_membre("exemple7 fois 2.png", code2.fct_test, z), z
    )
    u = code2.inverse_csr(M[0], M[1])
    IB = M[2]
    j = x // dx
    i = y // dy
    if (i, j) in IB:
        return u[IB.index((i, j))]
    else:
        return 0


def f(x, y, z=1):
    d = code2.f_d(z)
    return P(x, y, z) - np.cos(2 * np.pi / d * y) * np.exp(-complex(0, k * x))


verification(1)

L1_max, L2_max, H1_max = 0, 0, 0
INF=[]
L1=[]
L2=[]
H1=[]




# for r in range(3):
    
#     z = [0.5, 1, 2][r]
#     fichier = ["exemple7 divise 2.png", "exemple7.png", "exemple7 fois 2.png"][r]
#     matrice_reacteur = code2.Mapping.mapping_BC(code2.Mapping.color_to_flag(code2.Mapping.mapping(code2.Mapping.png_to_rgb_matrix(fichier))))
#     n,p=matrice_reacteur.shape
#     dx=code2.f_dx(z)
#     dy=code2.f_dy(z)
#     d=code2.f_d(z)
#     M=code2.matrice_to_csr(matrice_reacteur, chi, code2.second_membre(matrice_reacteur,code2.fct_test,z),z)
#     u=code2.inverse_csr(M[0],M[1])
#     IB=M[2]
#     x_range=[0,n*dx]
#     y_range=[0,p*dy]

#     data = np.zeros((p, p))


#     max_difference = -float("inf")  #norme infini
#     somme=0                         #norme L1 **2
#     somme_carres=0                  #norme L2 **2
#     somme_carres1=0                 #norme L2 de la dérivée **2
#     f_reference = []
#     df_ref=0
    
#     for x in np.arange(x_range[0], x_range[1], dx):
#         j=int(x//dx)
#         for y in np.arange(y_range[0], y_range[1], dy):
#             i=p-int(y//dy)-1
#             f_reference.append(abs(np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k * x))))
            
#             if (i,j) in IB:
#                 P=u[IB.index((i,j))]
#                 difference = P - np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k * x))
#                 somme+=abs(difference)
#                 somme_carres+=abs(P - np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k * x)))**2
#                 data[j,i]=difference
#             else:
#                 difference=0
#                 data[j,i]=difference
#             if j>0 and i>0:
#                     if (i,j+1) in IB and (i,j-1) in IB:
#                         derivee_x=(u[IB.index((i,j+1))]-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x+dx)))-u[IB.index((i,j-1))]+np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x-dx))))/(2*dx)
#                         df_x=(-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x+dx)))+np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x-dx))))/(2*dx)
#                     elif (i,j+1) in IB:
#                         derivee_x= (u[IB.index((i,j+1))]-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x+dx))))/(2*dx)
#                         df_x=(-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x+dx))))/(2*dx)
#                     elif (i,j-1) in IB:
#                         derivee_x= (-u[IB.index((i,j-1))]+np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x-dx))))/(2*dx)
#                         df_x=(np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x-dx))))/(2*dx)
#                     else:
#                         derivee_x=0
#                         df_x=0
#                     if (i+1,j) in IB and (i-1,j) in IB:
#                         derivee_y=(u[IB.index((i-1,j))]-np.cos(2 * np.pi / d * (y+dy)) * np.exp(-complex(0,k*x))-u[IB.index((i+1,j))]+np.cos(2 * np.pi / d * (y-dy)) * np.exp(-complex(0,k*x)))/(2*dy)
#                         df_y=(-np.cos(2 * np.pi / d * (y+dy)) * np.exp(-complex(0,k*x))+np.cos(2 * np.pi / d * (y-dy)) * np.exp(-complex(0,k*x)))/(2*dy)
#                     elif (i+1,j) in IB:
#                         derivee_y=(-u[IB.index((i+1,j))]+np.cos(2 * np.pi / d * (y-dy)) * np.exp(-complex(0,k*x)))/(2*dy)
#                         df_y=(np.cos(2 * np.pi / d * (y-dy)) * np.exp(-complex(0,k*x)))/(2*dy)
#                     elif (i-1,j) in IB:
#                         derivee_y=(-np.cos(2 * np.pi / d * (y+dy)) * np.exp(-complex(0,k*x))+u[IB.index((i-1,j))])/(2*dy)
#                         df_y=(-np.cos(2 * np.pi / d * (y+dy)) * np.exp(-complex(0,k*x)))/(2*dy)
#                     else:
#                         derivee_y=0
#                         df_y=0
#                     somme_carres1+=abs(derivee_x)**2+abs(derivee_y)**2
#                     df_ref+=abs(df_x)**2+abs(df_y)**2

#             if abs(difference) > max_difference:
#                 max_difference = abs(difference)
#     # print(somme)
#     # print((sum([el for el in f_reference])))
#     print(r)
#     INF.append((max_difference) / max(f_reference))
#     L1.append(somme / sum(f_reference))
#     L2.append((np.sqrt(somme_carres)) / np.sqrt(sum([el**2 for el in f_reference])))
#     H1.append((np.sqrt(somme_carres)+np.sqrt(somme_carres1)) / (np.sqrt(sum([el**2 for el in f_reference]))+np.sqrt(df_ref)))

# plt.scatter(['précision fois 1', "précision fois 2", "précision fois 4"] ,L1,color='red',label='norme L1')
# plt.scatter(['précision fois 1', "précision fois 2", "précision fois 4"] ,L2,color='green',label='norme L2')
# plt.scatter(['précision fois 1', "précision fois 2", "précision fois 4"] ,INF,color='black',label='norme infini')
# plt.scatter(['précision fois 1', "précision fois 2", "précision fois 4"], H1,color='blue',label='norme H1')
# plt.legend()
# plt.title('normes pour 3 précisions différentes')

# # # Créez une figure avec un seul sous-graphique
# # fig, ax = plt.subplots()

# # # Affichez les données avec imshow
# # im = ax.imshow(data, cmap='coolwarm')

# # # Ajoutez une barre de couleur
# # cbar = plt.colorbar(im)
# plt.show()
