
function j = BilinearQuadElementConductiveMatrix(h,k_darcy,x1,y1,x2,y2,x3,y3,x4,y4)

syms s t;
a = (y1*(s-1)+y2*(-1-s)+y3*(1+s)+y4*(1-s))/4;
b = (y1*(t-1)+y2*(1-t)+y3*(1+t)+y4*(-1-t))/4;
c = (x1*(t-1)+x2*(1-t)+x3*(1+t)+x4*(-1-t))/4;
d = (x1*(s-1)+x2*(-1-s)+x3*(1+s)+x4*(1-s))/4;

B1 = [a*(t-1)/4-b*(s-1)/4 ; c*(s-1)/4-d*(t-1)/4];

B2 = [a*(1-t)/4-b*(-1-s)/4; c*(-1-s)/4-d*(1-t)/4];

B3 = [a*(t+1)/4-b*(s+1)/4; c*(s+1)/4-d*(t+1)/4];

B4 = [a*(-1-t)/4-b*(1-s)/4 ; c*(1-s)/4-d*(-1-t)/4];

Bfirst = [B1 B2 B3 B4];
Jfirst = [0 1-t t-s s-1 ; t-1 0 s+1 -s-t ;
  s-t -s-1 0 t+1 ; 1-s s+t -t-1 0];
J = [x1 x2 x3 x4]*Jfirst*[y1 ; y2 ; y3 ; y4]/8;%此处是计算行列式J的值
B = Bfirst/J;

j0=J*k_darcy*transpose(B)*B;
j1=int(int(j0, t, -1, 1), s, -1, 1);
j=double(j1)*h;