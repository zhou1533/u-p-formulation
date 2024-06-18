clear;
clc;
syms s;
E=3e3;
NU=0.2;
h=1;%单元厚度，若h取1，就可以不考虑h
n1=0.333;%这里表示孔隙度
density_solid=0.306;%定义密度等于2kg/m3
density_water=0.2977;
k_darcy=0.004883;
g=10;
Kf=3.999e4;
Ks=1e25;
Qb=1/((n1/Kf+(1-n1)/Ks));%1/Qb=(n1/Kf+(a-n1)/Ks) 这里a为固相骨架可压缩系数，a=1-Kd/Ks,因Ks>>Kd所以，a=1.Kd,Ks,Kf分别表示固体骨架的体积模量、土颗粒的体积模量以及流体体积模量
k_dynamic=k_darcy/(density_water*g);%动态达西渗透系数k_dynamic=k_darcy/(desity_water*g)



for i=1:101%这个i表示的是土柱高度 
    i1=i-1;
  for  j=1:1%这个循环表示的无量纲时间的循环
   % tao=20+(j-1)*20;
   tao=60;
    
lamda=833.3;
alpha=1;
Mu=1250;
Vc=sqrt((lamda+2*Mu+alpha^2*Qb)/density_solid);
beta=density_water/density_solid;
kappa=Qb/(lamda+2*Mu+alpha^2*Qb);
kama=beta/n1;
b=1/(kappa-beta^2);
a=(kama-beta^2)/b;



temp1=alpha/((1-alpha*beta)*sqrt(a));

fun=@(s)(exp(-b*s/(2*a))).*besseli(0,b*sqrt(s.^2-a*i1.^2)/(2*a)).*heaviside(s-i1*sqrt(a));

W(i,j)=temp1*vpa(int(fun,s,0.0,tao));
    
fun1=@(s) heaviside(s-i1);%这里没有用f（tao）相当于f(tao)=1;

temp5=vpa(int(fun1,s,0.0,tao));

u(i,j)=-temp5-beta*W(i);
  end
end



figure (1)
plot(0:1:100,u(:,1),'g','linewidth', 1.5)
% hold on 
% plot(0:1:100,u(:,2),'b','linewidth', 1.5)
% hold on
% plot(0:1:100,u(:,3),'r','linewidth', 1.5)








