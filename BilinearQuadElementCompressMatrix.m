function y = BilinearQuadElementCompressMatrix(h,Qb,x1,y1,x2,y2,x3,y3,x4,y4)
%BilinearQuadElementMass  

syms s t;

N1=(1-s)*(1-t)/4;
N2=(1+s)*(1-t)/4;
N3=(1+s)*(1+t)/4;
N4=(1-s)*(1+t)/4;

N=[N1,N2,N3,N4];

Jfirst = [0 1-t t-s s-1 ; t-1 0 s+1 -s-t ;
  s-t -s-1 0 t+1 ; 1-s s+t -t-1 0];
J = [x1 x2 x3 x4]*Jfirst*[y1 ; y2 ; y3 ; y4]/8;%此处是计算行列式J的值


s0=J*(1/Qb)*transpose(N)*N;
s1=int(int(s0, t, -1, 1), s, -1, 1);
y= double(h*s1);
