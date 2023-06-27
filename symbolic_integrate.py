import sympy as sym
a0,a1,a2=sym.symbols('a0,a1,a2')
f=a0*a0+a1*a1+a2*a2+(1-a0-a1-a2)*(1-a0-a1-a2)
int0=sym.integrate(f,(a0,0,1))
print(int0)
int1=sym.integrate(int0,(a1,0,1))
int2=sym.integrate(int1,(a2,0,1))
print(int2)

