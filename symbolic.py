import sympy as sym
f00,f01,f02,f10,f11,f12,f20,f21,f22=sym.symbols('f00,f01,f02,f10,f11,f12,f20,f21,f22')
C,D=sym.symbols('C,D')
F=sym.Matrix([[f00,f01,f02],[f10,f11,f12],[f20,f21,f22]])
J=F.det()
traceFTF=(F.transpose()*F).trace()
psi=C*((J**(-sym.Rational(2,3)))*traceFTF-3)+D*((J-1)*(J-1))
print('dw(0)=', sym.ccode(sym.diff(psi,f00)),";")
print('dw(1)=', sym.ccode(sym.diff(psi,f10)),";")
print('dw(2)=', sym.ccode(sym.diff(psi,f20)),";")
print('dw(3)=', sym.ccode(sym.diff(psi,f01)),";")
print('dw(4)=', sym.ccode(sym.diff(psi,f11)),";")
print('dw(5)=', sym.ccode(sym.diff(psi,f21)),";")
print('dw(6)=', sym.ccode(sym.diff(psi,f02)),";")
print('dw(7)=', sym.ccode(sym.diff(psi,f12)),";")
print('dw(8)=', sym.ccode(sym.diff(psi,f22)),";")

df0=sym.diff(psi,f00)
df1=sym.diff(psi,f10)
df2=sym.diff(psi,f20)
df3=sym.diff(psi,f01)
df4=sym.diff(psi,f11)
df5=sym.diff(psi,f21)
df6=sym.diff(psi,f02)
df7=sym.diff(psi,f12)
df8=sym.diff(psi,f22)
print("Result for d2vdq2:")
df_list=[df0,df1,df2,df3,df4,df5,df6,df7,df8]
var_list=[f00,f10,f20,f01,f11,f21,f02,f12,f22]
for r in range(9):
    for c in range(9):
        print('ddw('+str(r)+','+str(c)+')=',sym.ccode(sym.diff(df_list[r],var_list[c])),';')


