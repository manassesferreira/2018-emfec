library(menina)
g=6
mu=10
r=0.2
eps=0.02
D=definaDispositivo(g,mu)

a=c()
for(d in seq(0,length(D$seg))){
    a=c(a,paste(D$seg[d],D$pos[d]))
}
fileConn<-file("in.txt")
writeLines(a, fileConn)
close(fileConn)

C=computeComponente(D,g,r,eps)
B=busqueBase(C,D,g,r,eps)
A=assenteAcesso(B,C,D,g,r,eps)


print(length(unique(C$comp)))
print(length(A$x))
for(a in seq(0,length(A$x))){
    print(paste(A$x[a],A$y[a]))
}
