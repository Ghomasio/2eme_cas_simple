
#libraries to import
from matplotlib import pyplot as plt
import numpy as np

#discrétisation

Longueur_x = 1
Nb_x = 30
dx = Longueur_x/(Nb_x)

Longueur_y = 1
Nb_y = 30
dy = Longueur_y/(Nb_y)

#ecriture de la matrice (matrice A)

mat = np.zeros((Nb_x*Nb_y, Nb_x*Nb_y))          #La fonction np.eye() de NumPy crée une matrice diagonale où les éléments de la diagonale sont égaux à 1 et tous les autres éléments sont égaux à zéro.
mat[:Nb_y + 1, :Nb_y + 1] = np.eye(Nb_y + 1)    #définit les conditions aux limites pour les bords gauche et droit du domaine. Elle affecte les Nb_y + 1 premières lignes et colonnes de mat en tant qu'identité, ce qui signifie que les valeurs sur ces bords sont fixées.
mat[-Nb_y:, -Nb_y:] = np.eye(Nb_y)              #fait de même pour les bords supérieur et inférieur du domaine


for i in range(1, Nb_x+1): #colonnes
    for j in range(1, Nb_y+1): #lignes

        k = j + (i-1)*Nb_y - 1

        #TEST CONDITIONS LIMITES
        if j == 1:                        # bord du bas
            mat[k, k] = -3/(2*dy)
            mat[k, k+1] = 4/(2*dy)
            mat[k, k+2] = -1 /(2*dy) 
        
        elif i == Nb_x :                        # bord droite
            mat[k, k] = 3/(2*dx)
            mat[k, k-Nb_x] = -4/(2*dx)
            mat[k, k-2*Nb_x] = 1 /(2*dx)  

        elif i == 1 or j == Nb_y:      #Si le point est sur un bord, il est déjà traité par les conditions aux limites, donc la valeur de la matrice à cet emplacement est fixée à 1
            mat[k, k] = 1           # bord gauche et bord du haut

        #CORPS DE LA MATRICE    
        else : 
            mat[k, k - Nb_x] = mat[k, k + Nb_x] = 1/(dx**2)     #On remplie les termes à gauche et à droite 
            mat[k, k - 1] = mat[k, k + 1] = 1/(dy**2)           #On remplie les termes en haut et en bas
            mat[k, k] = -2/(dx**2) - 2/(dy**2)                  #On remplie le terme centrale

#Conditions initiales(matrice B)

conditions_limites = np.zeros((Nb_x*Nb_y))

for i in range(1, Nb_x+1):
    conditions_limites[i*Nb_y-1] = np.sin(np.pi*i*dx)




#Résolution du sytème AX = B

resol_temp = np.linalg.inv(mat)@conditions_limites



#reconstruction la solution sous forme de grille bidimensionnelle

#parcourt chaque élément de la solution temporaire resol_temp et le place à 
#la position correspondante dans la grille bidimensionnelle sol. Les indices sont 
#inversés (-j-1, -i) car la solution est stockée dans l'ordre inverse par rapport 
#à la manière dont les boucles sont structurées. Cela permet d'obtenir une représentation 
#correcte de la solution dans la grille bidimensionnelle sol.

sol = np.zeros((Nb_x, Nb_y))

for i in range(Nb_y) :
    for j in range(Nb_x):
        sol[-j-1, i-1] = resol_temp[j + (i-1)*Nb_y]


#CHARGE HYDRAULIQUE :


h = sol


# Générer une grille pour positionner les vecteurs
x = np.linspace(0, Longueur_x, Nb_x)
y = np.linspace(Longueur_y, 0, Nb_y)
X, Y = np.meshgrid(x, y)



#on sait que v = -kgrad(h)

gradient_x = np.zeros((Nb_x, Nb_y))
gradient_y = np.zeros((Nb_x, Nb_y))

Vx = np.zeros((Nb_x, Nb_y))
Vy = np.zeros((Nb_x, Nb_y))

k = 10e-5                  #Pour notre sujet

def compute_gradient(matrix):
    gradient_x = np.zeros_like(matrix, dtype=float)
    gradient_y = np.zeros_like(matrix, dtype=float)

    # Calculer le gradient selon l'axe x
    gradient_x[:, 1:-1] = (matrix[:, 2:] - matrix[:, :-2]) / 2.0
    gradient_x[:, 0] = matrix[:, 1] - matrix[:, 0]
    gradient_x[:, -1] = matrix[:, -1] - matrix[:, -2]

    # Calculer le gradient selon l'axe y
    gradient_y[1:-1, :] = (matrix[2:, :] - matrix[:-2, :]) / 2.0
    gradient_y[0, :] = matrix[1, :] - matrix[0, :]
    gradient_y[-1, :] = matrix[-1, :] - matrix[-2, :]

    return gradient_x, gradient_y


gradient_x , gradient_y = compute_gradient(h)

Vx = -k * gradient_x
Vy = k * gradient_y

# Tracer la carte de charge hydraulique
plt.figure(figsize=(10, 8))
plt.contourf(X, Y, h, cmap='viridis', levels=25)
plt.colorbar(label='Charge hydraulique')

# Tracer les vecteurs de vitesse
plt.quiver(X, Y, Vx, Vy, scale=0.0001, color='white')

plt.xlabel('Position x')
plt.ylabel('Position y')
plt.title('Charge hydraulique et Vecteurs de vitesse')
plt.show()

