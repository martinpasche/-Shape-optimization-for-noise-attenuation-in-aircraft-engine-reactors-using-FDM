import math
import numpy as np
import sympy
import code2
import adjoint
import matplotlib.pyplot as plt

k0=code2.k0
Z0=code2.Z0
Z=code2.Z
alpha=k0*Z0/Z
M0=code2.M0

def energie(u,T):
    #T=code2.Mapping.mapping_BC(code2.Mapping.color_to_flag(code2.Mapping.mapping(code2.Mapping.png_to_rgb_matrix(image))))
    # u=code2.solution(T,g,chi,z)[0]
    u=u[0]
    somme=0
    n,p=T.shape[0],T.shape[1]
    for k in range(0,n*p):
        (i,j)=code2.Mapping.matricise(k,p)
        if T[i,j]['flag'] in ['EDP','frontiere']:
            somme+=abs(u[i,j])**2
    return somme

def derivee_energie(p,v,dp,dv,T,chi):
    p=p[0]
    #T=code2.Mapping.mapping_BC(code2.Mapping.color_to_flag(code2.Mapping.mapping(code2.Mapping.png_to_rgb_matrix(image))))
    n,P=T.shape[0],T.shape[1]
    somme=0
    for k in range(0,n*P):
        (i,j)=code2.Mapping.matricise(k,P)
        if 'BC' in T[i,j] and T[i,j]['BC']=='Neumann':
            somme+=(alpha.imag*(p[i,j].real*v[i,j].real-p[i,j].imag*v[i,j].imag)+2*alpha.real*(p[i,j].imag*v[i,j].real+p[i,j].real*v[i,j].imag)-2*alpha.real*M0*(dp[i,j].real*v[i,j].real-dp[i,j].imag*v[i,j].imag)/k0+2*alpha.imag*M0*(dp[i,j].real*v[i,j].imag+dp[i,j].imag*v[i,j].real)/k0+alpha.real*(M0/k0)**2*(dp[i,j].imag*dv[i,j].real-dp[i,j].real*dv[i,j].imag)+alpha.imag*(M0/k0)**2*(dp[i,j].real*dv[i,j].real+dp[i,j].imag*dv[i,j].imag))*chi(i,j)
    return somme

def dichotomie(n,p,T,chi0,beta,f,a,b,e):
    while b-a>e:
        m = (b+a)/2
        if f(n,p,T,chi0,beta,b)*f(n,p,T,chi0,beta,m)>0:
            b = m
        else:
            a = m
    return m

def equation(n,p,T,chi0,beta,l):
        integrale_chi=0
        for k in range(0,n*p):
            (i,j)=code2.Mapping.matricise(k,p)
            if 'BC' in T[i,j] and T[i,j]['BC']=='Neumann':
                integrale_chi+=max(0,min(chi0(i,j)+l,1))
        return integrale_chi-beta

def distance(mat1,mat2):
    return np.sum(np.abs(mat1 - mat2))

def f_adjoint(u,x,y,T,z):
    u=u[0]
    dy=adjoint.f_dy(z)
    dx=adjoint.f_dx(z)
    n,p=T.shape
    i=int(y//dy)
    j=int(x//dx)
    return (-2*u[i,j].real,2*u[i,j].imag)

def optimisation(chi00, zeta0, beta, delta,T,g,f,z=1):
    r = 0
    #T=code2.Mapping.mapping_BC(code2.Mapping.color_to_flag(code2.Mapping.mapping(code2.Mapping.png_to_rgb_matrix(image))))
    n,P=T.shape[0],T.shape[1]
    chi0_matrix = np.zeros((n,P))
    chi0=lambda i,j:chi00(T,i,j)
    for k in range(0,n*P):
        (i,j)=code2.Mapping.matricise(k,P)
        chi0_matrix[i,j]=chi0(i,j)
    solution = dichotomie(n,P,T,lambda x, y: chi0_matrix[x, y],beta,equation,0,1,delta)   #dichotomie pour trouver la solution l qui permet de définir la projection qui vérifie le volume de matériau donné
    #dp=code2.d_solution_dx(T,g,lambda x, y: chi0_matrix[x, y],z)
    
    chi1_matrix=np.zeros((n,P))
    p=code2.solution(T,g,lambda x, y: chi0_matrix[x, y],z)
    v=adjoint.solution(p,T,f,lambda x, y: chi0_matrix[x, y],z)
    dp=code2.d_solution_dx(T,g,lambda x, y: chi0_matrix[x, y],z)
    dv=adjoint.d_adjoint_dx(p,T,f,lambda x, y: chi0_matrix[x, y],z)
    for k in range(0,n*P):
        (i,j)=code2.Mapping.matricise(k,P)
        chi1_matrix[i,j]= max(0, min(chi0_matrix[i,j]-zeta0*derivee_energie(p,v,dp,dv,T,chi0)+solution, 1))
    
    while distance(chi0_matrix,chi1_matrix) >= delta:
        u0=code2.solution(T,g,lambda x, y: chi0_matrix[x, y],z)
        u1=code2.solution(T,g,lambda x, y: chi1_matrix[x, y],z)
        r = r+1
        if energie(u0,T)>energie(u1,T):
            zeta0=zeta0+0.001
                
        chi0_matrix=chi1_matrix
    
        solution = dichotomie(n,P,T,lambda x, y: chi0_matrix[x, y],beta,equation,0,1,delta)
        #dp=code2.d_solution_dx(T,g,lambda x, y: chi0_matrix[x, y],z)
        p=code2.solution(T,g,lambda x, y: chi0_matrix[x, y],z)
        v=adjoint.solution(p,T,f,lambda x, y: chi0_matrix[x, y],z)
        dp=code2.d_solution_dx(T,g,lambda x, y: chi0_matrix[x, y],z)
        dv=adjoint.d_adjoint_dx(p,T,f,lambda x, y: chi0_matrix[x, y],z)
        for k in range(0,n*P):
            (i,j)=code2.Mapping.matricise(k,P)
            chi1_matrix[i,j]= max(0, min(chi0_matrix[i,j]-zeta0*derivee_energie(p,v,dp,dv,T,lambda x, y: chi0_matrix[x, y])+solution, 1))
        
        else:
            k=0
            while energie(u0,T)<energie(u1,T):
                zeta0=zeta0/2**(k+1)
                chi0_matrix=chi1_matrix
        
                solution = dichotomie(n,P,T,lambda x, y: chi0_matrix[x, y],beta,equation,0,1,delta)
                #dp=code2.d_solution_dx(T,g,lambda x, y: chi0_matrix[x, y],z)
                p=code2.solution(T,g,lambda x, y: chi0_matrix[x, y],z)
                v=adjoint.solution(p,T,f,lambda x, y: chi0_matrix[x, y],z)
                dp=code2.d_solution_dx(T,g,lambda x, y: chi0_matrix[x, y],z)
                dv=adjoint.d_adjoint_dx(p,T,f,lambda x, y: chi0_matrix[x, y],z)
                for k in range(0,n*P):
                    (i,j)=code2.Mapping.matricise(k,P)
                    chi1_matrix[i,j]= max(0, min(chi0_matrix[i,j]-zeta0*derivee_energie(p,v,dp,dv,T,lambda x, y: chi0_matrix[x, y])+solution, 1))
                u0=code2.solution(T,g,lambda x, y: chi0_matrix[x, y],z)
                u1=code2.solution(T,g,lambda x, y: chi1_matrix[x, y],z)
                k+=1
    return chi1_matrix



chi=optimisation(code2.chi,0.1,0.5,0.1,code2.matrice_reacteur,code2.fct_test,f_adjoint,1)

min = np.min(chi)
max = np.max(chi)

fig, (ax1, ax2) = plt.subplots(1, 2)

# Afficher la partie réelle dans le premier sous-graphique
im1 = ax1.imshow(chi, cmap="coolwarm", interpolation="none")
ax1.set_title("chi optimal non projeté")
im1.set_clim(vmin=min, vmax=max)  # Ajuster l'échelle de couleur
fig.colorbar(im1, ax=ax1)
plt.show()
