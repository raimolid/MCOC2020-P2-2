  
import numpy as np

g = 9.81 #kg*m/s^2


class Barra(object):

    """Constructor para una barra"""
    def __init__(self, ni, nj, R, t, E, ρ, σy):
        super(Barra, self).__init__()
        self.ni = ni
        self.nj = nj
        self.R = R
        self.t = t
        self.E = E
        self.ρ = ρ
        self.σy = σy

    def obtener_conectividad(self):
        return [self.ni, self.nj]

    def calcular_area(self):
        A = np.pi*(self.R**2) - np.pi*((self.R-self.t)**2)
        return A

    def calcular_largo(self, reticulado):
        """Devuelve el largo de la barra. 
        ret: instancia de objeto tipo reticulado
        """
        xi = reticulado.obtener_coordenada_nodal(self.ni)
        xj = reticulado.obtener_coordenada_nodal(self.nj)
        dij = xi-xj
        return np.sqrt(np.dot(dij,dij))

    def calcular_peso(self, reticulado):
        """Devuelve el largo de la barra. 
        ret: instancia de objeto tipo reticulado
        """
        L = self.calcular_largo(reticulado)
        A = self.calcular_area()
        return self.ρ * A * L * g
    
    def obtener_rigidez(self, ret):
        L = self.calcular_largo(ret)
        A = self.calcular_area()
        k = self.E * A / L
        Tϴ = np.matrix([-np.cos(60), -np.sin(60), np.cos(60), np.sin(60)])
        ke = Tϴ.T @ Tϴ * k
        return np.array(ke)
        
    def obtener_vector_de_cargas(self, ret):
        W = self.calcular_peso(ret)
        v = np.array([[0, -1, 0, -1]])
        fe = (v.T) * W/2
        return fe
        
    def obtener_fuerza(self, ret):
        L = self.calcular_largo(ret)
        A = self.calcular_area()
        
        ni = self.ni
        nj = self.nj
        ue = np.array([ret.u[ni*2], ret.u[((2*ni)+1)], ret.u[2*nj], ret.u[((2*nj)+1)]]) 
        
        
        se = ((A * self.E)/L) * np.array([-np.cos(60), -np.sin(60), np.cos(60), np.sin(60)].T) * ue 
        return se
        
        

       
        
       