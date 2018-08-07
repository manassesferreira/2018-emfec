#include "common.h"
#include "existeConectividade.h"
#include "segmentosDoCruzamento.h"
#include "weighted_quick_union.h"

#include <Rcpp.h>
using namespace Rcpp;

NumericVector unaEdescubra(int Vertices, NumericVector u, NumericVector v) {
  NumericVector out(Vertices);
  int Edges = u.size();
  if(Edges>0){
    UF *uf;
    uf = new weightedQuickUnion(Vertices);
    int node1, node2;
    //cout<<Edges<<" "<<v.size()<<endl;
//    cout << " union "  << endl;
    for(int i = 0; i < Edges; ++i) {
      //cout << i << ": " << u[i] << "," << v[i] << endl;
      node1=u[i]; node2=v[i];
      uf->_union(node1, node2);
    }
//    cout << " find "  << endl;
    for (int i = 0; i < Vertices; i++) {
      //cout << i << ": " << uf->find(i) << endl;
      out[i] = uf->find(i);
/*
      for (int j = i+1; j < Vertices; j++) {
        if (uf->connected(i, j)){
          cout << "(" << i << "," << j << ") "<< uf->connected(i, j) << endl;
        }
      }
*/
    }
  }
  return out ;
}

// [[Rcpp::export]]
List avalieAcesso(List A, List D, int g, double r, double eps) {
 int Dispositivos = as<NumericVector>(D[0]).size();
 int Segmentos = 2 * g * ( g - 1 );
 int mu=Dispositivos/Segmentos;
 int s;

 NumericVector seg = as<NumericVector>(D[0]);
 NumericVector pos = as<NumericVector>(D[1]);

 int Acessos = as<NumericVector>(A[0]).size();
 //cout << Acessos << endl;

 NumericVector Aseg = as<NumericVector>(A[0]);
 NumericVector Apos = as<NumericVector>(A[1]);

 NumericVector _x(Dispositivos+Acessos);
 NumericVector _y(Dispositivos+Acessos);


 typedef struct {
   std::list<int> lista;
 } TipoArrayList;
 TipoArrayList *dispositivosNoSegmento;
 dispositivosNoSegmento=new TipoArrayList[Segmentos];

 for (int i = 0; i < Dispositivos; ++i) {
   s = (int)seg[i];
   //std::cout << i << " " << s << " " << pos[i] << std::endl;
   dispositivosNoSegmento[s].lista.push_back(i);
   if( s%(2*g-1) > (g-2) ){ //segmento vertical
     //std::cout << "\t" << (double)(s%(2*g-1)-(g-1)) << " " << pos[i] << std::endl;
     _x[i] =  (double)(s%(2*g-1)-(g-1)) ;
     _y[i] =  pos[i] + s/(2*g-1) ;
   }else{ //segmento horizontal
     //std::cout << "\t" << pos[i] << " " << (double)(floor(s/(2*g-1))) << std::endl;
     _x[i] =  pos[i] + s%(2*g-1) ;
     _y[i] =  (double)(floor(s/(2*g-1))) ;
   }
 }

 for (int i = Dispositivos; i < Dispositivos+Acessos; ++i) {
   s = (int)Aseg[i-Dispositivos];
   //std::cout << i << " " << s << " " << Apos[i-Dispositivos] << std::endl;
   dispositivosNoSegmento[s].lista.push_back(i);
   if( s%(2*g-1) > (g-2) ){ //segmento vertical
     _x[i] =  (double)(s%(2*g-1)-(g-1)) ;
     _y[i] =  Apos[i-Dispositivos] + s/(2*g-1) ;
   }else{ //segmento horizontal
     _x[i] =  Apos[i-Dispositivos] + s%(2*g-1) ;
     _y[i] =  (double)(floor(s/(2*g-1))) ;
   }
 }

 //cout << " _u _v "  << endl;

 int NedgesMAX=Segmentos*mu+6*g*g + Acessos*4  + Acessos*(Acessos-1)/2; // 2*g*(g-1)*mu+4*g*mu + 6*g*g
 int NsitesCRUZ=4*g*g; // 4*g*g

 NumericVector _u(NedgesMAX);
 NumericVector _v(NedgesMAX);

 NumericVector _minDoSeg(NsitesCRUZ);
 double minSeg;
 NumericVector _maxDoSeg(NsitesCRUZ);
 double maxSeg;

 int edge_counter=0;

 //nos segmentos
 for (int i =0; i < 2*(g-1)*g; ++i){

   double minSeg=1.1;
   int dmin=-1;
   double maxSeg=-0.1;
   int dmax=-1;


   //ordenando a lista pela posicao do dispositivo, do menor para o maior
   std::vector<std::pair<double,int>> auxiliar; 
   std::list<int>::const_iterator d1;
   for (d1 = dispositivosNoSegmento[i].lista.begin();
     d1 != dispositivosNoSegmento[i].lista.end(); ++d1) {
     std::pair<double, int> p = {pos[*d1], *d1}; 
     auxiliar.push_back(p);
   }
   std::sort(auxiliar.begin(), auxiliar.end());

   std::list<int>::iterator di;
   di = dispositivosNoSegmento[i].lista.begin();
   for ( std::pair<double,int> & element : auxiliar) {
//       cout << i << " " << element.first << " " << element.second << endl; 

//       cout << "antes " << *di << endl;
       *di=element.second;
//       cout << "depois " << *di << endl;

       ++di;
   }

   std::list<int>::const_iterator aux,d2;//mas mantenha o respeito
   for (d1 = dispositivosNoSegmento[i].lista.begin();
     d1 != dispositivosNoSegmento[i].lista.end(); ++d1) {
     aux=d1;
     d2=++aux; 
     if( d2 != dispositivosNoSegmento[i].lista.end() ) {
       if ( fabs(pos[*d1] - pos[*d2]) <= r ){
         _u[edge_counter]=*d1;
         _v[edge_counter]=*d2;
         edge_counter=edge_counter+1;
       }
     }
     if ( pos[*d1] > maxSeg ){
        maxSeg=pos[*d1];
        dmax=*d1;
     }
     if ( pos[*d1] < minSeg ){
        minSeg=pos[*d1];
        dmin=*d1;
     }
   }

   //estou assumindo que ha pelo menos um dispositivo no segmento
   _maxDoSeg[i] = dmax;
   _minDoSeg[i] = dmin;

 }


 //nos cruzamentos
 for (int i = 0; i < g*g; ++i) {
   //std::cout << i << " : ";
   int borda; //1-OS,2-S,3-LS,4-O,5-L,6-NO,7-N,8-NL
   int N, L, O, S;
   int xCruz = i%g;
   int yCruz = floor(i/g);

   //cout << "\t segmentosDoCruzamento "  << endl;
   segmentosDoCruzamento(g,xCruz,yCruz,&N,&L,&O,&S,&borda);

   //novamente assumindo que teremos pelo menos um dispositivo por segmento
   std::list<int> dispositivosNoCruzamento;
   dispositivosNoCruzamento.push_back(_minDoSeg[N]);
   dispositivosNoCruzamento.push_back(_minDoSeg[L]);
   dispositivosNoCruzamento.push_back(_maxDoSeg[O]);
   dispositivosNoCruzamento.push_back(_maxDoSeg[S]);

   //std::copy(std::begin(dispositivosNoCruzamento), std::end(dispositivosNoCruzamento),
   //        std::ostream_iterator<int>(std::cout, " "));
   //std::cout << std::endl;
   //std::cout << std::endl;

   std::list<int>::const_iterator d1;
   double d1x, d1y;
   //cout << "\t d1 d2 "  << endl;
   std::list<int>::const_iterator d2,aux;
   double d2x, d2y;
   for (d1 = dispositivosNoCruzamento.begin();
     d1 != dispositivosNoCruzamento.end(); ++d1) {
     aux = d1;
     for (d2 = ++aux; d2 != dispositivosNoCruzamento.end(); ++d2) {
       //std::cout << *d1 <<  " " << *d2 << std::endl;

       d1x = _x[*d1]; d1y = _y[*d1]; 
       //std::cout << "\t" << d1x <<  " " << d1y << std::endl;

       d2x = _x[*d2]; d2y = _y[*d2];
       //std::cout << "\t" << d2x <<  " " << d2y << std::endl;

       if ( existeConectividade(g, borda, r, eps, d1x, d1y, d2x, d2y) ) {
         //std::cout << "\t" << *d1 <<  " " << *d2 << std::endl;
         _u[edge_counter]=*d1;
         _v[edge_counter]=*d2;
         edge_counter=edge_counter+1;
       }
     }
   }
   //std::cout << std::endl;
 }
 
 for(int i = Dispositivos; i < Dispositivos+Acessos; ++i){
   for(int j = i+1; j < Dispositivos+Acessos; ++j){
     _u[edge_counter]=i;
     _v[edge_counter]=j;
     edge_counter=edge_counter+1;
   }
 }

// cout << " _C "  << endl;
  NumericVector _C = unaEdescubra(Dispositivos+Acessos, _u, _v); //using union-find method

 delete [] dispositivosNoSegmento;
 return Rcpp::List::create(Rcpp::Named("x") = _x,
                           Rcpp::Named("y") = _y,
                           Rcpp::Named("u") = _u,
                           Rcpp::Named("v") = _v,
                           Rcpp::Named("comp") = _C);
}
