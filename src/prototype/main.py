import random
import itertools
from itertools import groupby

#==============================================================================#
def segmentosDoCruzamento(g, xCruz, yCruz):
    LCruz = (2*g-1)*yCruz + xCruz
    NCruz = LCruz + (g-1)
    OCruz = LCruz - 1
    SCruz = LCruz - g

    if xCruz == 0:
        OCruz=LCruz + (g-2);
        if yCruz == 0: #*borda=1;
            SCruz=(2*g-1)*(g-1)-g + xCruz
        elif yCruz == (g-1): #*borda=6;
            NCruz=xCruz + (g-1)
        #else: #*borda=4;
            #nothing
    elif xCruz == (g-1):
        LCruz=LCruz - (g-1)
        if yCruz == 0: #*borda=3;
            SCruz=(2*g-1)*(g-1)-g + xCruz
        elif yCruz==(g-1):    #*borda=8;
            NCruz=xCruz + (g-1)
        #else: #*borda=5;
            #else
    else:
        if yCruz==0: #*borda=2;
            SCruz=(2*g-1)*(g-1)-g + xCruz
        elif yCruz==(g-1): #*borda=7;
            NCruz=xCruz + (g-1)
        #else #*borda=0;
            #matters
    return [NCruz, LCruz, OCruz, SCruz]

def dispositivosNoCruzamento(g, xCruz, yCruz):
    NLOS = segmentosDoCruzamento(g, xCruz, yCruz)
    sN = NLOS[0]
    sL = NLOS[1]
    sO = NLOS[2]
    sS = NLOS[3]

    dN = sN * mu
    dL = sL * mu
    dO = (sO + 1 ) * mu - 1
    dS = (sS + 1 ) * mu - 1

    return [dN, dL, dO, dS]

def computeComponentes(n_D, Edges):
    parent = [None] * n_D
    size = [None] * n_D
    for i in range(0,n_D):
        parent[i] = i
        size[i] = 1

    def find(p):
        root = p
        while root != parent[root]:
            root = parent[root]
        while p != root:
            newp = parent[p]
            parent[p] = root
            p = newp
        return root

    def union(p, q):
        rootP = find(p)
        rootQ = find(q)
        if rootP == rootQ:
            return
        if size[rootP] < size[rootQ]:
            parent[rootP] = rootQ
            size[rootQ] += size[rootP]
        else:
            parent[rootQ] = rootP
            size[rootP] += size[rootQ]

    for (u,v) in Edges:
        union(u, v)

    for (u,v) in Edges:
        union(u, v)

    C = list()
    labels=list(set(parent))
    for l in labels:
        C.append([i for i, e in enumerate(parent) if e == l])

    return C, parent

#==============================================================================#
g=10
a=1
mu=10

#==============================================================================#
n_s = 2 * g * (g - 1)
n_D = n_s * mu

Ds = [None] * n_D
for s in range(n_s):
    for i in range(mu):
        d = s * mu + i
        Ds[d] = s

#Dq = [None] * n_D
#for d in range(n_D):
#        Dq[d] = random.uniform(0, a)

Dq = [None] * n_D
for s in range(n_s):
    ru = [None] * mu
    for i in range(mu):
        ru[i] = random.uniform(0, a)
    ru = sorted(ru)
    for i in range(mu):
        d = s * mu + i
        Dq[d] = ru[i]

#==============================================================================#
T = 2 * g - 1
Dx = [None] * n_D
Dy = [None] * n_D
for d in range(n_D):
    s = Ds[d]
    q = Dq[d]
    if s % T > g - 2 :
        #print(s, "ver")
        Dx[d] = s % T - (g - 1)
        Dy[d] = q + s / T * a
    else:
        #print(s, "hor")
        Dx[d] = q + s % T * a
        Dy[d] = s / T * a

#==============================================================================#
L = (g - 1) * a

R = 0.25 * a
Dr = [R] * n_D
#for d in range(n_D):
#        Dr[d] = random.uniform(0, a)

epsilon = 0.02 * a / 2
delta = 2 * epsilon


n_A_max = n_s * (mu - 1) + g * g * 6
Arestas = [None] * n_A_max
Edges = list()
n_a = 0

for s in range(n_s):
    #print(s)
    #dmin = s * mu + 0
    #dmax = s * mu + mu - 1
    #print(dmin, Dq[dmin], dmax, Dq[dmax])

    for i in range(mu-1):
        d1 = s * mu + i
        d2 = d1 + 1
        Rmin=Dr[d1]
        if Dr[d2] < Rmin:
            Rmin=Dr[d2]
            #print(s, d1, d2, Dq[d2] - Dq[d1], Rmin)
        if Dq[d2] - Dq[d1] < Rmin:
            #print "seg", d1, d2
            Arestas[n_a] = (d1, d2)
            n_a = n_a +1
            Edges.append((d1, d2))

for xCruz in range(g):
    for yCruz in range(g):
        dNLOS = dispositivosNoCruzamento(g, xCruz, yCruz)
        #print xCruz, yCruz, dNLOS

        for (d1, d2) in itertools.combinations(dNLOS, 2):
            dx = abs(Dx[d2] - Dx[d1])
            dx = min(dx, L - dx)
            dy = abs(Dy[d2] - Dy[d1])
            dy = min(dy, L - dy)

            Rmin = Dr[d1]
            if Dr[d2] < Rmin:
                Rmin = Dr[d2]

            if dx == 0:
                if dy < Rmin:
                    #print "cruz1", d1, d2
                    Arestas[n_a] = (d1, d2)
                    n_a = n_a +1
                    Edges.append((d1, d2))
            elif dy == 0:
                if dx < Rmin:
                    #print "cruz2", d1, d2
                    Arestas[n_a] = (d1, d2)
                    n_a = n_a +1
                    Edges.append((d1, d2))
            else:
                if dx * dx + dy * dy < Rmin * Rmin:
                    if (
                        dx < epsilon or dy < epsilon or ( dx < delta and dy < delta ) or
                        ( dx < delta and abs(dy * (1 - epsilon / dx)) < epsilon ) or
                        ( dy < delta and abs(dx * (1 - epsilon / dy)) < epsilon )
                        ):
                        #print "cruz3", d1, d2
                        Arestas[n_a] = (d1, d2)
                        n_a = n_a +1
                        Edges.append((d1, d2))

#==============================================================================#
Componentes, Dc = computeComponentes(n_D, Edges)

#==============================================================================#
Bases = list()
Bx = list()
By = list()
Bd = list()

labels=list(set(Dc))
for l in labels:
    d=Dc.index(l)
    base=[l]
    if base not in Bases:
        Bases.append(base)
        Bx.append(Dx[d])
        By.append(Dy[d])
        Bd.append([d])

for s in range(n_s):
    #print(s)
    for i in range(mu-1):
        d1 = s * mu + i
        d2 = d1 + 1
        Rmin = Dr[d1]
        if Dr[d2] < Rmin:
            Rmin = Dr[d2]
        if Dc[d2] != Dc[d1]:
            if Dq[d2] - Dq[d1] < 2 * Rmin:
                #print "seg base", d1, d2, (Dq[d2] + Dq[d1]) / 2, Dc[d1], Dc[d2]
                base=[Dc[d1],Dc[d2]]
                if base not in Bases:
                    Bases.append(base)
                    #print base
                    Bs = s
                    Bq = (Dq[d2] + Dq[d1]) / 2
                    if Bs % T > g - 2 :
                        x = Bs % T - (g - 1)
                        y = Bq + Bs / T * a
                    else:
                        x = Bq + Bs % T * a
                        y = Bs / T * a
                    Bx.append(x)
                    By.append(y)
                    Bd.append([d1, d2])

for xCruz in range(g):
    for yCruz in range(g):
        dNLOS = dispositivosNoCruzamento(g, xCruz, yCruz)
        #print xCruz, yCruz, dNLOS

        B=list()
        Ds=list()
        for (d1, d2) in itertools.combinations(dNLOS, 2):
            if Dc[d2] != Dc[d1]:
                dx = abs(Dx[d2] - Dx[d1])
                dx = min(dx, L - dx)
                dy = abs(Dy[d2] - Dy[d1])
                dy = min(dy, L - dy)

                Rmin = Dr[d1]
                if Dr[d2] < Rmin:
                    Rmin = Dr[d2]

                if dx == 0:
                    if dy < 2 * Rmin:
                        #print "cruz1 base", d1, d2
                        B.append(Dc[d1])
                        B.append(Dc[d2])
                        Ds.append([d1, d2])
                elif dy == 0:
                    if dx < 2 * Rmin:
                        #print "cruz2 base", d1, d2
                        B.append(Dc[d1])
                        B.append(Dc[d2])
                        Ds.append([d1, d2])
                else:
                    if Dx[d2] == round(Dx[d2]):
                        Lx = abs(Dx[d1] - round(Dx[d1]))
                        Ny = abs(Dy[d2] - round(Dy[d2]))
                    else:
                        Lx = abs(Dx[d2] - round(Dx[d2]))
                        Ny = abs(Dy[d1] - round(Dy[d1]))

                    if Lx > Ny:
                        Bq = max(Lx - Rmin, 0)
                        dB = Bq * Bq + Ny * Ny
                        dBx = Bq
                        dBy = Ny
                    else:
                        Bq = max(Ny - Rmin, 0)
                        dB = Lx * Lx + Bq * Bq
                        dBx = Lx
                        dBy = Bq

                    if dB < Rmin * Rmin:
                        if (
                            dBx < epsilon or dBy < epsilon or ( dBx < delta and dBy < delta ) or
                            ( dBx < delta and abs(dBy * (1 - epsilon / dBx)) < epsilon ) or
                            ( dBy < delta and abs(dBx * (1 - epsilon / dBy)) < epsilon )
                            ):
                            #print "cruz3 base", d1, d2
                            B.append(Dc[d1])
                            B.append(Dc[d2])
                            Ds.append([d1, d2])
        #print B
        Bf = [len(list(group)) for key, group in groupby(sorted(B))]
        #print Bf

#        if len(Bf) > 3:
#            if min(Bf) > 0:
#                #print("quatro-tres-dois")
#                base=frozenset(B)
#                Bases.add(base)
#                #print base
#                Bx.append(xCruz)
#                By.append(yCruz)
#        elif len(Bf) > 2:
#            if min(Bf) > 0:
#                #print("tres-dois")
#                base=frozenset(B)
#                Bases.add(base)
#                #print base
#                Bx.append(xCruz)
#                By.append(yCruz)
        if len(Bf) > 1:
            if min(Bf) > 0:
                #print("dois")
                base=list(set(B))
                if base not in Bases:
                    Bases.append(base)
                    #print base
                    Bx.append(xCruz)
                    By.append(yCruz)
                    Bd.append(Ds)

#==============================================================================#
#solucao trivial
def trivial(Componentes):
    n_B = len(Componentes)
    map_MC2B = [None] * n_B
    for i in range(0,n_B):
        map_MC2B[i] = i
    return n_B, map_MC2B

#solucao monte carlo
def monteCarlo(escolhidos, universo):
    map_MC2B = random.sample(range(0,universo), k=escolhidos)
    return escolhidos, map_MC2B


#==============================================================================#
#n_B, map_MC2B = trivial(Componentes)
#n_B, map_MC2B = monteCarlo(1,len(Bx))

def passoMC(Nb, Bases, Edges, n_D):
    nB, map_MC2B = monteCarlo(Nb,len(Bases))

    n_MC=n_D+n_B
    Edges_MC = list()
    for e in Edges:
        Edges_MC.append(e)

    for b1 in range(n_D,n_MC-1):
        b2 = b1 + 1
        Edges_MC.append((b1,b2))

    for b1 in range(n_D,n_MC):
        for d in Bd[map_MC2B[b1-n_D]]:
            if isinstance(d, int):
                Edges_MC.append((b1, d))
            else:
                if len(d)>1:
                    Edges_MC.append((b1, d[0]))
                    Edges_MC.append((b1, d[1]))
                else:
                    Edges_MC.append((b1, d[0]))

    Componentes_MC, Dc_MC = computeComponentes(n_MC, Edges_MC)
    return Componentes_MC, map_MC2B

#==============================================================================#
def selecione(quantos, possiveis):
    return random.sample(possiveis, k=quantos)

def avalie(n_D, Edges, n_B, map_MC2B):
    n_MC = n_D + n_B
    Edges_MC = list()
    for e in Edges:
        Edges_MC.append(e)

    for b1 in range(n_D, n_MC-1):
        b2 = b1 + 1
        Edges_MC.append((b1,b2))

    for b1 in range(n_D, n_MC):
        for d in Bd[map_MC2B[b1-n_D]]:
            if isinstance(d, int):
                Edges_MC.append((b1, d))
            else:
                if len(d)>1:
                    Edges_MC.append((b1, d[0]))
                    Edges_MC.append((b1, d[1]))
                else:
                    Edges_MC.append((b1, d[0]))

    Componentes_MC, Dc_MC = computeComponentes(n_MC, Edges_MC)
    return Componentes_MC

#==============================================================================#
#print len(Componentes), Componentes
#print "---"
passos_MC=1000 #len(Componentes)

NumeroDeBases = len(Componentes)
##print "N_pa\t", NumeroDeBases
escolhidos_min = list()
N_min = len(Componentes)
n_B, map_MC2B = trivial(Componentes)
universo = range(0, len(Bases))
##print "|U|\t", len(universo)
escolhidos = list()
N_comp_esc = len(Componentes)
Found = False
Nb = 1
while Nb < NumeroDeBases+1 and not Found:
    print "\t", escolhidos_min
    possiveis = list(set(universo)-set(escolhidos_min)-set(escolhidos))
    candidatos = list(selecione(1, possiveis))
    for e in escolhidos_min:
        candidatos.append(e)

    for step in range(0,passos_MC):
#        print step, candidatos
        Componentes_MC = avalie(n_D, Edges, Nb, candidatos)
        N_comp_MC = len(Componentes_MC)
#        print "\t", N_previous, N_comp_MC


        if N_comp_esc > N_comp_MC:
            N_comp_esc = N_comp_MC
            escolhidos = list()
            for c in candidatos:
                escolhidos.append(c)

        else:
            if random.uniform(0, 1) < float(N_comp_esc) / N_comp_MC:
                N_comp_esc = N_comp_MC
                escolhidos = list()
                for c in candidatos:
                    escolhidos.append(c)

        if N_comp_esc < N_min:
            N_min = N_comp_esc
            escolhidos_min = list()
            for c in candidatos:
                escolhidos_min.append(c)
 
            #Componentes_MC = avalie(n_D, Edges, Nb, escolhidos_min)
            #print len(Componentes_MC), Componentes_MC
 
            if N_min == 1:
                Found = True

        possiveis = list(set(universo)-set(candidatos))
        candidatos = list(selecione(1, possiveis))
        if len(escolhidos_min) < Nb:
            for e in escolhidos_min:
                candidatos.append(e)
        else:
            for e in list(selecione(Nb-1, escolhidos_min)):
                candidatos.append(e)



    Nb = Nb + 1

print ""
if Found:
    print "Found!"
    print "\t", escolhidos_min
    print "\t", Nb-1, "ao inves de", NumeroDeBases, "de um universo de", len(universo)

else:
    print "Not found..."
    print "\t", N_min, " ao inves de 1"
