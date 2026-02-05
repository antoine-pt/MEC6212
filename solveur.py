import numpy as np
import sympy as sp


class Parametres():
    """Classe contenant les parametres du problème."""
    def __init__(self, nPts, R, D, S, Ce):
        self.nPts = nPts
        self.R = R
        self.D = D
        self.S = S
        self.Ce = Ce
        
    @property
    def dr(self):
        return self.R / (self.nPts - 1)
    
    @property
    def pos(self):
        return np.linspace(0, self.R, self.nPts)
    

def normeL1(ana, sim):
    """ Calcule la norme L1 entre une solution analytique et une solution
    numérique.

    Args:
        ana (np.array): solution analytique
        sim (np.array): solution numérique

    Returns:
        float: norme L1 entre les deux solutions
    """
    return np.sum(np.abs(ana - sim)) / np.size

def normeL2(ana, sim):
    """ Calcule la norme L2 entre une solution analytique et une solution
    numérique.

    Args:
        ana (np.array): solution analytique
        sim (np.array): solution numérique

    Returns:
        float: norme L2 entre les deux solutions
    """
    return np.sqrt(np.sum((ana - sim)**2) / ana.size)

def normeLinf(ana, sim):
    """ Calcule la norme Linf entre une solution analytique et une solution
    numérique.

    Args:
        ana (np.array): solution analytique
        sim (np.array): solution numérique

    Returns:
        float: norme Linf entre les deux solutions
    """
    return np.max(np.abs(ana - sim))

def solveur(params):

    A = np.zeros((nPts,nPts))
    b = np.zeros(nPts)
    pos = np.linspace(0,R,nPts)

    for i in range(nPts):
        ri = pos[i]
        if i == 0:
            A[i,i] = -1
            A[i,i+1] = 1

            b[i] = 0

        elif i == nPts-1:
            A[i,i] = 1
            b[i] = params.Ce

        else:
            A[i,i-1] = 1 / params.dr**2
            A[i,i] = -1 / (params.dr*ri)- 2 / params.dr**2
            A[i,i+1] = 1 / (params.dr*ri)+ 1 / params.dr**2

            b[i] = params.S/ params.D

    return np.linalg.solve(A,b)

def analytique(params):
    """Fonction renvoyant le vecteur solution analytique au probleme en fonction des parametres."""
    ## Solution analytique
    r = sp.Symbol('r')
    func = params.S / (4 * params.D) * ( r**2 - params.R**2) + params.Ce 
    anal = sp.lambdify(r, func, 'numpy')

    return anal(params.pos)

if __name__ == "__main__":
    nPts = 5
    R = 0.5
    D = 1e-10
    S = 2e-8
    Ce = 20

    params = Parametres(nPts=nPts,R=R,D=D,S=S,Ce=Ce)
    sim = solveur(params)
    ana = analytique(params)

    print(params.dr)
    print(params.pos)
    params.nPts = 20

    print(params.dr)
    print(params.pos)

    

