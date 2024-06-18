%这个结果是准确的，准确性见Simon论文图4
%论文中用图的图例是直接在图上添加编辑的
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




  for  j=1:81 %这个j表示的是时间
  fprintf('j %f***循环次数***\r\n',j);      
    
  tao=(j-1)*2.5;
    i=0;
    
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

digits(4);
fun=@(s)((k_darcy/Vc).*exp(-b*s/(2*a))).*besseli(0,b*sqrt(s.^2-a*i.^2)/(2*a)).*heaviside(s-i*sqrt(a));
%fun=@(s)(exp(-b*s/(2*a))).*besseli(0,b*sqrt(s.^2-a*i.^2)/(2*a)).*heaviside(s-i*sqrt(a));

W(j)=single(temp1*vpa(int(fun,s,0.0,tao)));
    
fun1=@(s) (k_darcy/Vc)*heaviside(s-i);
digits(4);

%fun1=@(s) -1*heaviside(s-i);
temp5=single(vpa(int(fun1,s,0.0,tao)));

u(j)=-temp5-beta*W(j);
u(j)=u(j)*Vc/k_darcy;
 end

u'


figure (1)
plot(0:2.5:200,u,'--r','linewidth', 2.5)
xlabel('\tau=t/ρk','fontsize',15);
ylabel('$${\hat{u}}=uV_c/k\sigma_0$$','interpreter','latex','fontsize',15,'fontweight','bold')
%('\hat{u}=uV_c/k\sigma_0','interpreter','latex'); 

 