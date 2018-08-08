resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}

library("menina")
library("RColorBrewer")


lambdas=seq(10,10)
L=20
r=1/3
delta=1/20

passosTOTAL=1000
intervalo=passosTOTAL/10
passos=seq(1,passosTOTAL)
instancias=seq(1,1)


for( lambda in lambdas ){
  for( instancia in instancias ){
    N_seg=2*L*(L-1)
    N_dis=N_seg*lambda
 
    D=definaDispositivo(L,lambda)
    C=computeComponente(D,L,r,delta)
    N_comp=length(unique(C$comp))
 
    N_pa_MAX=N_dis
    if( N_comp < N_pa_MAX ){
      N_pa_MAX=N_comp
    }
 
    B=busqueBase(C,D,L,r,delta)
    A=assenteAcesso(B,C,D,L,r,delta)
    N_aprox=length(A$x)

    npas=seq(1,N_aprox)
 
    N_previous=N_comp;
    N_min=N_comp;
    for(npa in npas){
#      N_C=c()
#      step=c()
      for(passo in passos){

        seg=floor(runif(npa, 0, N_seg-1))
        pos=runif(npa, 0, 1)
        A_MC=data.frame( seg=seg, pos=pos )
        A_MCx=c()
        A_MCy=c()
        for(i in seq(1,npa)){
          s = A_MC$seg[i];
          if( s%%(2*L-1) > (L-2) ){ #segmento vertical
            A_MCx = c(A_MCx,  s%%(2*L-1)-(L-1) );
            A_MCy = c(A_MCy, A_MC$pos[i] + floor(s/(2*L-1)) );
          }else{ #segmento horizontal
            A_MCx = c(A_MCx, A_MC$pos[i] + s%%(2*L-1) );
            A_MCy = c(A_MCy, floor(s/(2*L-1)) );
          }
        }

        C_MC=avalieAcesso(A_MC,D,L,r,delta)
        N_comp_MC=length(unique(C_MC$comp))
 
        if ( N_previous > N_comp_MC) {
          N_previous = N_comp_MC ;
        }else{
          if(runif(1,0,1) < N_previous/N_comp_MC ){
            N_previous = N_comp_MC;
          }
        }
#        if ( passo %% intervalo == 1 ) {
#          N_C=c(N_C,N_previous)
#          step=c(step,passo)
#        }

        if(N_previous < N_min){
#          N_C=c(N_C,N_previous)
          step=c(step,passo)
          N_min = N_previous

          par(mar=c(0,0,0,0)+0.1)
          plot( -10,-10, xlim=c(0,L-1), ylim=c(0,L-1), xaxt='n', yaxt='n', ann=FALSE)
  
          labelsUNIQ=unique(C_MC$comp)
          colorsUNIQ=rainbow(length(labelsUNIQ))
          colorID=match(C_MC$comp,labelsUNIQ)
          colors=colorsUNIQ[colorID]
  
          text(C_MC$x,C_MC$y, labels=colorID, cex= 0.48, col=colors)
          points(A_MCx,A_MCy,pch=2,cex=2)

          points(A$x,A$y,pch=2,col=3,cex=2)
        }

      }
      print(paste("npa",npa,"N_min",N_min))
#      par(resetPar())
#      plot( step, N_C, type='l',ylim=c(0,max(N_C)),main=paste("Instancia",instancia,"; N_pa =",npa))
#      abline(h=N_aprox,col=3)


#      print("aprox")
#      Aseg=c()
#      Apos=c()
#      for(i in seq(1,N_aprox)){
#         Ax=A$x[i]
#         Ay=A$y[i]
         #print(paste(Ax,Ay))
#         if(Ax == floor(Ax)){
#            if(Ay == floor(Ay)){
               #print("vertical e horizontal") #inserindo como vertical...
#            }else{
               #print("vertical")
#            }
#            As=floor(Ay)*(2*L-1)+(L-1)+Ax
#            Ap=Ay-floor(Ay)
#         }else{
#             if(Ay == floor(Ay)){
               #print("horizontal")
#               As=Ay*(2*L-1)+floor(Ax)
#               Ap=Ax-floor(Ax)
#            }else{
               #print("wtf") #shit never happens
#            }        
#         }
         #print(paste(As,Ap))
#         Aseg=c(Aseg,As)
#         Apos=c(Apos,Ap)
#      }
#      A_aprox=data.frame(seg=Aseg,pos=Apos)
#      C_aprox=avalieAcesso(A_aprox,D,L,r,delta)
#      N_comp_aprox=length(unique(C_aprox$comp))
#      abline(h=N_comp_aprox,col=3,lty=3)

    }
  }
}
