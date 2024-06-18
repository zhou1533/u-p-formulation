function y = BilinearQuadAssemble_Q(Q,q,i,j,m,n)
%BilinearQuadAssemble      This function assembles the element
%                          
Q(2*i-1,i) = Q(2*i-1,i) + q(1,1);
Q(2*i-1,j) = Q(2*i-1,j) + q(1,2);
Q(2*i-1,m) = Q(2*i-1,m) + q(1,3);
Q(2*i-1,n) = Q(2*i-1,n) + q(1,4);

Q(2*i,i) = Q(2*i,i) + q(2,1);
Q(2*i,j) = Q(2*i,j) + q(2,2);
Q(2*i,m) = Q(2*i,m) + q(2,3);
Q(2*i,n) = Q(2*i,n) + q(2,4);

Q(2*j-1,i) = Q(2*j-1,i) + q(3,1);
Q(2*j-1,j) = Q(2*j-1,j) + q(3,2);
Q(2*j-1,m) = Q(2*j-1,m) + q(3,3);
Q(2*j-1,n) = Q(2*j-1,n) + q(3,4);

Q(2*j,i) = Q(2*j,i) + q(4,1);
Q(2*j,j) = Q(2*j,j) + q(4,2);
Q(2*j,m) = Q(2*j,m) + q(4,3);
Q(2*j,n) = Q(2*j,n) + q(4,4);

Q(2*m-1,i) = Q(2*m-1,i) + q(5,1);
Q(2*m-1,j) = Q(2*m-1,j) + q(5,2);
Q(2*m-1,m) = Q(2*m-1,m) + q(5,3);
Q(2*m-1,n) = Q(2*m-1,n) + q(5,4);

Q(2*m,i) = Q(2*m,i) + q(6,1);
Q(2*m,j) = Q(2*m,j) + q(6,2);
Q(2*m,m) = Q(2*m,m) + q(6,3);
Q(2*m,n) = Q(2*m,n) + q(6,4);

Q(2*n-1,i) = Q(2*n-1,i) + q(7,1);
Q(2*n-1,j) = Q(2*n-1,j) + q(7,2);
Q(2*n-1,m) = Q(2*n-1,m) + q(7,3);
Q(2*n-1,n) = Q(2*n-1,n) + q(7,4);

Q(2*n,i) = Q(2*n,i) + q(8,1);
Q(2*n,j) = Q(2*n,j) + q(8,2);
Q(2*n,m) = Q(2*n,m) + q(8,3);
Q(2*n,n) = Q(2*n,n) + q(8,4);

y = Q;