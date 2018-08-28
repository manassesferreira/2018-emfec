from menina import definaDispositivo
from menina import computeComponente
from menina import busqueBase
from menina import assenteAcesso

g=6
a=1
mu=10
r=0.2
R=r*a
n_s=2*g*(g-1)
n_D=n_s*mu
Dr=[R]*n_D
eps=0.02

passos_MC=100

Ds, Dq = definaDispositivo(g,a,mu)
C, Dc = computeComponente(r,eps,Dr,n_D,n_s,Ds,Dq,g,a,mu)
Bases, Bx, By, Bd = busqueBase(C,Dc,r,eps,Dr,n_D,n_s,Ds,Dq,g,a,mu)
Ax, Ay = assenteAcesso(passos_MC,Bases,Bx,By,Bd,C,Dc,r,eps,Dr,n_D,n_s,Ds,Dq,g,a,mu)

print len(set(Dc))
print len(Ax)
for a in range(len(Ax)):
    print Ax[a],Ay[a]
