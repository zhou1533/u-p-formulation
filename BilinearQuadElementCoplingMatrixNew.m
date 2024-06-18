function q = BilinearQuadElementCoplingMatrixNew(h,x1,y1,x2,y2,x3,y3,x4,y4)
%���Ӻ�����Ŀ���Ǽ�����Ͼ���Q�е�Ԫ���Ӿ���

syms s t;

N1=(1-s)*(1-t)/4;
N2=(1+s)*(1-t)/4;
N3=(1+s)*(1+t)/4;
N4=(1-s)*(1+t)/4;


N=[N1,0,N2,0,N3,0,N4,0;0,N1,0,N2,0,N3,0,N4];

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
J = [x1 x2 x3 x4]*Jfirst*[y1 ; y2 ; y3 ; y4]/8;%�˴��Ǽ�������ʽJ��ֵ
B = Bfirst/J;


q0=J*transpose(N)*B;%�ر�ע��Ҫ����J
q1=int(int(q0, t, -1, 1), s, -1, 1);
q=single(q1)*h;
