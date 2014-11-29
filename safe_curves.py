from edwards_curve import EdwardsCurve
from sage.rings.finite_rings.constructor import GF
from exceptions import TypeError

class SafeCurve(EdwardsCurve):
    def __init__(self, p, s, base_point_x, base_point_y, base_point_l):
        self._p = p
        self._K = GF(self._p)
        
        self._s = self._K(s)
        self._c = 2 / (self._s**2)
        self._r = self._c + (1/self._c)
        self._d = - ( (self._c + 1) ** 2 ) / ( (self._c - 1) ** 2 )
        
        self._base_point_x = base_point_x
        self._base_point_y = base_point_y
        self._base_point_l = base_point_l
        
        EdwardsCurve.__init__(self, self._K, self._d)
        
        self._base_point = self((self._base_point_x,self._base_point_y))
     
    def base_point(self):
        return self._base_point

def Curve1174():
    p = 2**251 - 9
    s = 1806494121122717992522804053500797229648438766985538871240722010849934886421
    base_x = 1582619097725911541954547006453739763381091388846394833492296309729998839514
    base_y = 3037538013604154504764115728651437646519513534305223422754827055689195992590
    base_l =  904625697166532776746648320380374280092339035279495474023489261773642975601
    return SafeCurve(p, s, base_x, base_y, base_l)    
        
        #=======================================================================
        # self._p = 2**251 - 9
        # self._K = GF(self._p)
        # 
        # self._s = self._K(1806494121122717992522804053500797229648438766985538871240722010849934886421)
        # self._d = self._K(-1174)
        # self._base_x = self._K(1582619097725911541954547006453739763381091388846394833492296309729998839514)
        # self._base_y = self._K(3037538013604154504764115728651437646519513534305223422754827055689195992590)
        # self._l =  904625697166532776746648320380374280092339035279495474023489261773642975601
        # 
        # 
        # EdwardsCurve.__init__(self, self._K, self._d)
        # 
        # self._base_point = self((self._base_x,self._base_y))
        #=======================================================================
        
class Curve25519(EdwardsCurve):
    def __init__(self):
        self.__p = 2**255 - 19
        self.__d = 486662
        EdwardsCurve.__init__(self,GF(self.__p),self.__d)

class CurveE222(EdwardsCurve):
    def __init__(self):
        self.__p = 2**222 - 117
        self.__d = 160102
        EdwardsCurve.__init__(self,GF(self.__p),self.__d)
        
class CurveE382(EdwardsCurve):
    def __init__(self):
        self.__p = 2**382 - 105
        self.__d = -67254
        EdwardsCurve.__init__(self,GF(self.__p),self.__d)
        
class Curve41417(EdwardsCurve):
    def __init__(self):
        self.__p = 2**251 - 9
        self.__d = -1174
        EdwardsCurve.__init__(self,GF(self.__p),self.__d)
        
class CurveEd448_Goldilocks(EdwardsCurve):
    def __init__(self):
        self.__p = 2**251 - 9
        self.__d = -1174
        EdwardsCurve.__init__(self,GF(self.__p),self.__d)
        
class CurveE521(EdwardsCurve):
    def __init__(self):
        self.__p = 2**251 - 9
        self.__d = -1174
        EdwardsCurve.__init__(self,GF(self.__p),self.__d)
        
