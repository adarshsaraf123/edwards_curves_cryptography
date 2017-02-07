from sage.crypto.cryptosystem import PublicKeyCryptosystem
from exceptions import TypeError, AttributeError
from sage.rings.integer import Integer
from sage.rings.arith import is_square
from sage.schemes.elliptic_curves.all import is_EllipticCurve
#from sage.schemes.elliptic_curves.edwards_curve import is_EdwardsCurve #to be added when you want to add this to the Sage library later
from sage.functions.other import sqrt
from sage.schemes.elliptic_curves import safe_curves

class EllipticCurveCryptosystem(PublicKeyCryptosystem):
    def __init__(self, E = safe_curves.Curve1174(), base_point = safe_curves.Curve1174().base_point()):
        if not (is_EllipticCurve(E) or is_EdwardsCurve(E)):
            raise TypeError("Elliptic Curve Cryptosystem is only over Elliptic curves or Edwards curves")
        self.__E = E
        self.__base_point = self._curve(base_point)
    
    def curve(self):
        return self.__E
    
    def base_point(self):
        return self.__base_point

class ECDH(EllipticCurveCryptosystem):
    def __init__(self, E = safe_curves.Curve1174(), base_point = safe_curves.Curve1174().base_point()):
        EllipticCurveCryptosystem.__init__(self, E, base_point)
   
    def decrypt(self, ciphertext_pair, B_private_key):
        first, second = ciphertext_pair
        return (second - (first * B_private_key) )
    def encrypt(self, plaintext, A_private_key, B_public_key):
        r"""
        For convenience of understanding we state the encryption process to be:
            Alice wants to encrypt the 'plaintext' using
                * Alice's private key: 'A_private_key'; it is an integer
                * Bob's public key: 'B_public_key'; it is a point on the elliptic curve on which the system is based
        Returns a tuple (A_private_key*curve.base_point(), self.encode(plaintext) + A_private_key*B_public_key)
        """
        #P_m = self.encode(plaintext) # plain message point
        P_m = plaintext
        try:
            P = self.curve()((B_public_key[0],B_public_key[1]))
        except TypeError:
            raise ValueError("Bob's public key not on the curve of the cryptosystem")
        
        try:
            k = Integer(A_private_key)
        except ValueError:
            raise ValueError("Alice's private key not a valid integer")
        
        try:
            G = self.base_point()
        except AttributeError:
            raise AttributeError("No base point defined for the curve")
        
        return (k*G, P_m + k*P) 
    
    def encode(self, message):
        """To be implemented for encoding a given bit-string as en elliptic curve point on the used curve for the cryptosystem
        """
        return self.curve()(message)
    
class Elligator():
    """To implement the Elligator map as presented in Bernstein, Lange (2011)
    """
    def __init__(self, E = safe_curves.Curve1174()):
        if not is_EdwardsCurve(E):
            raise AttributeError("Elligator is defined over Edwards curve only")
        self.__E = E
        
        self._s = self.get_s()
        if not self._s:
            raise ValueError("Elligator map is not defined for this Edwards Curve as -d is not a square, where d is the Edwards curve parameter")
               
        self._c = 2 / (self._s**2)
        self._r = self._c + (1/self._c)
        self._d = - ( (self._c + 1) ** 2 ) / ( (self._c - 1) ** 2 )
    
    def get_s(self):
        """Returns a single value of s as required for the Elligator map while there may even be two of them, in which case the returned values is the smaller of the two
        """ 
        try:
            return self._s
        except AttributeError:
            pass
        s = []
        d = self.__E.get_d()
        if is_square(-d):
            c = ( 2*(d-1) + 4*sqrt(-d) ) / ( 2*(d+1) )
            s_2 = 2/c
            if is_square(s_2):
                s.append(sqrt(s_2))
            
            c = ( 2*(d-1) - 4*sqrt(-d) ) / ( 2*(d+1) )
            s_2 = 2/c
            if is_square(s_2):
                s.append(sqrt(s_2))
        if len(s) == 0:
            return None
        elif len(s) == 1:
            return s[0]
        else:
            if s[0] < s[1]:
                return s[0]
            else: 
                return s[1]    
    
    def curve(self):
        return self.__E
    
    def chi(self, a):
        K = self.curve().base_ring()
        a = K(a)
        q = K.order() 
        return a ** ( (q-1) / 2 )
    
    def encode(self, t):
        K = self.curve().base_ring()
        t = K(t)
        q = K.order()
        s = self._s
        c = self._c
        r = self._r
        u = (1 - t)/(1 + t)
        v = u**5 + ( (r**2 - 2) * (u**3) ) + u
        X = self.chi(v) * u
        Y = ( (self.chi(v) * v) ** ((q+1)/4) ) * self.chi(v) * self.chi(u**2 + (1/c**2))
        x = (c-1) * s * X * (1+X) / Y
        y = (r*X - (1 + X)**2 ) / (r*X + (1 + X)**2 ) 
        return self.curve()( [x,y] )
 
    def is_image_point(self, point):
        s = self._s
        c = self._c
        r = self._r
        x = point[0]
        y = point[1]
        eta = (y-1) / (2 * (y+1))
        if y + 1 == 0:
            return False
        if not is_square( (1 + eta * r)**2 - 1):
            return False
        if eta*r == -2:
            temp = 2*s*(c - 1)*self.chi(c)/r
            if not x == temp:
                return False
        return True   
    
    def decode(self, point):
        point = self.curve()(point)
        K = self.curve().base_ring()
        q = K.order()
        s = self._s
        c = self._c
        r = self._r
        x = point[0]
        y = point[1]
        eta = (y-1) / (2 * (y+1))
        if self.is_image_point(point):
            X_bar = -(1 + eta*r) + ((1 + eta*r)**2 -1) ** ( (q+1)/4 )
            z = self.chi( (c - 1) * s * X_bar * (1 + X_bar) * x * (X_bar**2 + 1/c**2))
            u_bar = z * X_bar
            t_bar = (1 - u_bar) / (1 + u_bar)
        if t_bar > ((q-1)/2):
            t = -t_bar
        else:
            t = t_bar
        return t
    
    
    
    
    
    
    