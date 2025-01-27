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


n = 120
p = 30
L = 105  # longueur exemple 7 en pixel
h = 15  # hauteur exemple 7 en pixel
k = code2.k
M0 = code2.M0
x0 = 8
x1 = 112
y0 = 8
y1 = 22


def verification(z=1):
    dx = code2.f_dx(z)
    dy = code2.f_dy(z)
    d = code2.f_d(z)
    # Créez des valeurs pour x et y
    x = np.linspace(0, n, n)
    y = np.linspace(0, p, p)
    X, Y = np.meshgrid(x, y)
    print(X.shape)
    # Calculez votre fonction Z(x, y) - vous devez remplacer cette fonction par la vôtre
    # Exemple : une fonction simple
    Z = np.zeros((p, n), dtype=complex)
    condition = (X >= 8) & (X <= (112)) & (Y >= 8) & (Y <= 22)
    for i in range(0, p):
        for j in range(0, n):
            if condition[i, j] == True:
                Z[i,j] = (X[i,j]-x0)**2*(X[i,j]-x1)*(Y[i,j]-y0)**2*(Y[i,j]-y1)**2
                # Z[i,j] =-1/2*((Y[i,j]-y0)**2*(Y[i,j]-y1)**2*(k*k*((X[i,j]-x0)**2*(X[i,j]-x1))+2*(1-M0)*((X[i,j]-x1)+2*(X[i,j]-x0)))+\
                #     2*(X[i,j]-x0)**2*(X[i,j]-x1)*((Y[i,j]-y1)**2+(Y[i,j]-y0)**2+4*(Y[i,j]-y0)*(Y[i,j]-y1))\
                #     +complex(0,-1)*2*k*M0*(Y[i,j]-y0)**2*(Y[i,j]-y1)**2*(2*(X[i,j]-x0)*(X[i,j]-x1)+(X[i,j]-x0)**2))
                # np.cos(2 * np.pi / d * Y[i, j] * dy) * np.exp(
                #     -complex(0, k * X[i, j] * dx)
                # )
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
    ax1.set_title("Partie Réelle test")
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
    ax2.set_title("Partie Imaginaire test")
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
        matrice_reacteur, chi, code2.second_membre("exemple7.png", code2.fct_test, z), z
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

# INF=[]
# L1=[]
# L2=[]
# H1=[]
# for z in range(1,11):
#     dx=code2.f_dx(z)
#     dy=code2.f_dy(z)
#     d=code2.f_d(z)
#     M=code2.matrice_to_csr(matrice_reacteur, chi, code2.second_membre("exemple7.png",code2.fct_test,z),z)
#     u=code2.inverse_csr(M[0],M[1])
#     IB=M[2]
#     x_range=[0,n*dx]
#     y_range=[0,p*dy]


#     max_difference = -float("inf")  #norme infini
#     somme=0                         #norme L1 **2
#     somme_carres=0                  #norme L2 **2
#     somme_carres1=0                 #norme L2 de la dérivée **2
#     for x in np.arange(x_range[0], x_range[1], dx):
#         i=x//dx
#         for y in np.arange(y_range[0], y_range[1], dy):
#             j=y//dy
#             if (i,j) in IB:
#                 P=u[IB.index((i,j))]
#                 difference = P-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k * x))
#                 somme+=abs(P-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k * x)))
#                 somme_carres+=abs(P-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k * x)))**2
#             else:
#                 difference=-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k * x))
#                 somme+=abs(-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k * x)))
#                 somme_carres+=abs(-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k * x)))**2
#             if j>0 and i>0:
#                     if (i,j+1) in IB and (i,j-1) in IB:
#                         derivee_x=(u[IB.index((i,j+1))]-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x+dx)))-u[IB.index((i,j-1))]+np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x-dx))))/(2*dx)
#                     elif (i,j+1) in IB:
#                         derivee_x= (u[IB.index((i,j+1))]-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x+dx)))+np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x-dx))))/(2*dx)
#                     elif (i,j-1) in IB:
#                         derivee_x= (-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x+dx)))-u[IB.index((i,j-1))]+np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x-dx))))/(2*dx)
#                     else:
#                         derivee_x=(-np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x+dx)))+np.cos(2 * np.pi / d * y) * np.exp(-complex(0,k*(x-dx))))/(2*dx)
#                     if (i+1,j) in IB and (i-1,j) in IB:
#                         derivee_y=(u[IB.index((i+1,j))]-np.cos(2 * np.pi / d * (y+dy)) * np.exp(-complex(0,k*x))-u[IB.index((i-1,j))]+np.cos(2 * np.pi / d * (y-dy)) * np.exp(-complex(0,k*x)))/(2*dy)
#                     elif (i+1,j) in IB:
#                         derivee_y=(-np.cos(2 * np.pi / d * (y+dy)) * np.exp(-complex(0,k*x))-u[IB.index((i+1,j))]+np.cos(2 * np.pi / d * (y-dy)) * np.exp(-complex(0,k*x)))/(2*dy)
#                     elif (i-1,j) in IB:
#                         derivee_y=(-np.cos(2 * np.pi / d * (y+dy)) * np.exp(-complex(0,k*x))-u[IB.index((i-1,j))]+np.cos(2 * np.pi / d * (y-dy)) * np.exp(-complex(0,k*x)))/(2*dy)
#                     else:
#                         derivee_y=(-np.cos(2 * np.pi / d * (y+dy)) * np.exp(-complex(0,k*x))+np.cos(2 * np.pi / d * (y-dy)) * np.exp(-complex(0,k*x)))/(2*dy)
#                     somme_carres1+=abs(derivee_x)**2+abs(derivee_y)**2

#             if abs(difference) > max_difference:
#                 max_difference = abs(difference)
#     INF.append(max_difference)
#     L1.append(np.sqrt(somme))
#     L2.append(np.sqrt(somme_carres))
#     H1.append(np.sqrt(somme_carres)+np.sqrt(somme_carres1))
# plt.scatter(range(1,11),L1,color='red',label='norme L1')
# plt.scatter(range(1,11),L2,color='green',label='norme L2')
# plt.scatter(range(1,11),INF,color='black',label='norme infini')
# #plt.scatter(range(1,11),H1,color='blue',label='norme H1')
# plt.legend()
# plt.title('normes en fonction de z=0.01/dx avec dx en m')
# plt.show()
