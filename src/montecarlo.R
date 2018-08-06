library(menina)

lambdas=seq(10,10)
L=20
r=1/3
delta=1/20

passos=100

for( lambda in lambdas ){
  N_seg=L*(L-1)
  N_dis=N_seg*lambda

  D=definaDispositivo(L,lambda)
  C=computeComponente(D,L,r,delta)
  N_comp=length(unique(C$comp))

  N_pa_MAX=N_dis
  if( N_comp < N_pa_MAX ){
    N_pa_MAX=N_comp
  }

  for(npa in seq(1,N_pa_MAX)){
    for(passo in seq(1,passos)){
      print(paste(lambda,npa,passo))   

      seg=floor(runif(npa, 0, N_seg))
      pos=runif(npa, 0, 1)
      A_MC=data.frame( seg=seg, pos=pos )
      print(A_MC)
      #D_MC=data.frame( seg=c(D$seg,seg), pos=c(D$pos,pos) )

      C_MC=avalieAcesso(A_MC,D,L,r,delta)
      N_comp_MC=length(unique(C_MC$comp))

      print(N_comp_MC)

      print("--")
    }
  }

  #B=busqueBase(C,D,L,r,delta)
  #A=assenteAcesso(B,C,D,L,r,delta)
}
