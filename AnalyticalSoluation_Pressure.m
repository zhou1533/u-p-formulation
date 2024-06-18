
%求解结果应该如文章"AN ANALYTICAL SOLUTION FOR THE TRANSIENT RESPONSE OF SATURATED
%POROUS ELASTIC SOLIDS"中图2(b)。

%图例是直接在图上标注的

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




  for  j=1:81%这个j表示的是时间
  
    fprintf('j %f***循环次数***\r\n',j);  
%    if(j==1)
%       tao=0.0669;
%    else
%     tao=(j-1)*2.0078;
%    end
 tao=(j-1)*2.5;
 i=5;
    
    
lamda=833.3;
alpha=1;
Mu=1250.0;
Vc=sqrt((lamda+2*Mu+alpha^2*Qb)/density_solid);
beta=density_water/density_solid;
kappa=Qb/(lamda+2*Mu+alpha^2*Qb);
kama=beta/n1;
b=1/(kappa-beta^2);
a=(kama-beta^2)/b;



temp1=-alpha/(1-alpha*beta);


temp2=(k_darcy/Vc)*exp(-b*i/(2*sqrt(a))).*heaviside(tao-i*sqrt(a));
temp3=b*i/sqrt(a);
digits(6);%用于定义vpa的精度，4表示小数点后四位

fun=@(s)(k_darcy/Vc).*exp(-b*s/(2*a)).*besseli(1,b*sqrt(s.^2-a*i.^2)/(2*a)).*heaviside(s-i*sqrt(a))./sqrt(s.^2-a*i.^2);

W_e(j)=-1*(temp1*temp2+temp1*temp3*vpa(int(fun,s,0.0,tao)));%-1是考虑受力f(t)=-1;
    

u_e(j)=-(-1)*(k_darcy/Vc)*heaviside(tao-i)-beta*W_e(j);

p(j)=Vc*(beta*u_e(j)+kappa*W_e(j))/k_darcy;

 if(p(j)>=1)%%这里是擅自做的处理，因为计算结果有些不收敛
     p(j)=1.0;
 end


 end

p

figure (2)
%plot(0.0069:2.0078:200,p,'--r','linewidth', 1.5)
plot(0:2.5:200,p,'--r','linewidth', 2.5)
xlabel('\tau=t/ρk','fontsize',15);
 ylabel('$${\hat{p}}=p/\sigma_0$$','interpreter','latex','fontsize',15,'fontweight','bold')
gca=legend({'Analytical Solution'},'FontSize',14,'Location','SouthEast');
po=get(gca,'Position');%获取legend的位置
set(gca,'Position',[po(1)-0.05,po(2)+0.2,po(3),po(4)]);

hold on
