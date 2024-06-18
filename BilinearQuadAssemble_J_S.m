function y = BilinearQuadAssemble_J_S(J,jmatrix,i,j,m,n)
%BilinearQuadAssemble      This function assembles the element
%                          
J(i,i) = J(i,i) + jmatrix(1,1);
J(i,j) = J(i,j) + jmatrix(1,2);
J(i,m) = J(i,m) + jmatrix(1,3);
J(i,n) = J(i,n) + jmatrix(1,4);


J(j,i) = J(j,i) + jmatrix(2,1);
J(j,j) = J(j,j) + jmatrix(2,2);
J(j,m) = J(j,m) + jmatrix(2,3);
J(j,n) = J(j,n) + jmatrix(2,4);

J(m,i) = J(m,i) + jmatrix(3,1);
J(m,j) = J(m,j) + jmatrix(3,2);
J(m,m) = J(m,m) + jmatrix(3,3);
J(m,n) = J(m,n) + jmatrix(3,4);

J(n,i) = J(n,i) + jmatrix(4,1);
J(n,j) = J(n,j) + jmatrix(4,2);
J(n,m) = J(n,m) + jmatrix(4,3);
J(n,n) = J(n,n) + jmatrix(4,4);

 y= J;