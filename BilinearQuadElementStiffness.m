function w = BilinearQuadElementStiffness(E,NU,h,x1,y1,x2,y2,x3,y3,x4,y4,p)
%BilinearQuadElementStiffness       This function returns the element
%                                   stiffness matrix for a bilinear
%                                   quadrilateral element with modulus
%                                   of elasticity E, Poisson’s ratio
%                                   NU, thickness h, coordinates of
%                                   node 1 (x1,y1), coordinates
%node 3 (x3,y3), and coordinates of node 4 (x4,y4). Use p = 1 for casesof plane stress, and p = 2 for
%cases of plane strain.The size of the element stiffness matrix is 8 x 8.
syms s t;
a = (y1*(s-1)+y2*(-1-s)+y3*(1+s)+y4*(1-s))/4;
b = (y1*(t-1)+y2*(1-t)+y3*(1+t)+y4*(-1-t))/4;
c = (x1*(t-1)+x2*(1-t)+x3*(1+t)+x4*(-1-t))/4;
d = (x1*(s-1)+x2*(-1-s)+x3*(1+s)+x4*(1-s))/4;

B1 = [a*(t-1)/4-b*(s-1)/4 0 ; 0 c*(s-1)/4-d*(t-1)/4 ;
  c*(s-1)/4-d*(t-1)/4 a*(t-1)/4-b*(s-1)/4];

B2 = [a*(1-t)/4-b*(-1-s)/4 0 ; 0 c*(-1-s)/4-d*(1-t)/4;
  c*(-1-s)/4-d*(1-t)/4 a*(1-t)/4-b*(-1-s)/4];

B3 = [a*(t+1)/4-b*(s+1)/4 0 ; 0 c*(s+1)/4-d*(t+1)/4 ;
  c*(s+1)/4-d*(t+1)/4 a*(t+1)/4-b*(s+1)/4];

B4 = [a*(-1-t)/4-b*(1-s)/4 0 ; 0 c*(1-s)/4-d*(-1-t)/4 ;
  c*(1-s)/4-d*(-1-t)/4 a*(-1-t)/4-b*(1-s)/4];

Bfirst = [B1 B2 B3 B4];
Jfirst = [0 1-t t-s s-1 ; t-1 0 s+1 -s-t ;
  s-t -s-1 0 t+1 ; 1-s s+t -t-1 0];
J = [x1 x2 x3 x4]*Jfirst*[y1 ; y2 ; y3 ; y4]/8;%此处是计算行列式J的值
B = Bfirst/J;
if p == 1
  D = (E/(1-NU*NU))*[1, NU, 0 ; NU, 1, 0 ; 0, 0, (1-NU)/2];
elseif p == 2
  D = (E/(1+NU)/(1-2*NU))*[1-NU, NU, 0 ; NU, 1-NU, 0 ; 0, 0,(1-2*NU)/2];
end
BD = J*transpose(B)*D*B;%特别需要注意要乘J
r = int(int(BD, t, -1, 1), s, -1, 1);
z = h*r;
w = single(z);