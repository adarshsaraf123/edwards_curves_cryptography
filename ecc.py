from sage.crypto.cryptosystem import PublicKeyCryptosystem
import safe_curves
from exceptions import TypeError, AttributeError
from sage.rings.integer import Integer
from sage.rings.arith import is_square

class EllipticCurveCryptosystem(PublicKeyCryptosystem):
    def __init__(self, curve = safe_curves.Curve1174(), base_point = safe_curves.Curve1174().base_point()):
        self._curve = curve
        self._base_point = self._curve(base_point)
    
    def curve(self):
        return self._curve
    
    def base_point(self):
        return self._base_point

class ECDH(EllipticCurveCryptosystem):
    def __init__(self, curve = safe_curves.Curve1174(), base_point = safe_curves.Curve1174().base_point()):
        EllipticCurveCryptosystem.__init__(self, curve, base_point)
   
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
            raise TypeError("Bob's public key not on the curve of the cryptosystem")
        try:
            k = Integer(A_private_key)
        except TypeError:
            raise TypeError("Alice's private key not a valid integer")
        try:
            G = self.base_point()
        except AttributeError:
            raise AttributeError("No base point defined for the curve")
        
        return (k*G, P_m + k*P) 
    
    def encode(self, message):
        return self.curve()(message)
    
class Elligator():
    def __init__(self, curve = safe_curves.Curve1174()):
        self._curve = curve
    
    def curve(self):
        return self._curve
    
    def chi(self, a):
        K = self.curve().base_ring()
        a = K(a)
        q = K.order() 
        return a ** ( (q-1) / 2 )
    
    def encode(self, t):
        K = self.curve().base_ring()
        t = K(t)
        q = K.order()
        s = self.curve()._s
        c = self.curve()._c
        r = self.curve()._r
        u = (1 - t)/(1 + t)
        v = u**5 + ( (r**2 - 2) * (u**3) ) + u
        X = self.chi(v) * u
        Y = ( (self.chi(v) * v) ** ((q+1)/4) ) * self.chi(v) * self.chi(u**2 + (1/c**2))
        x = (c-1) * s * X * (1+X) / Y
        y = (r*X - (1 + X)**2 ) / (r*X + (1 + X)**2 ) 
        return self.curve()( [x,y] )
 
    def is_image_point(self, point):
        s = self.curve()._s
        c = self.curve()._c
        r = self.curve()._r
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
        s = self.curve()._s
        c = self.curve()._c
        r = self.curve()._r
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
    
    
    
    
    
    
    