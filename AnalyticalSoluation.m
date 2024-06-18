clear;
clc;
syms s;
E=3e3;
NU=0.2;
h=1;%��Ԫ��ȣ���hȡ1���Ϳ��Բ�����h
n1=0.333;%�����ʾ��϶��
density_solid=0.306;%�����ܶȵ���2kg/m3
density_water=0.2977;
k_darcy=0.004883;
g=10;
Kf=3.999e4;
Ks=1e25;
Qb=1/((n1/Kf+(1-n1)/Ks));%1/Qb=(n1/Kf+(a-n1)/Ks) ����aΪ����Ǽܿ�ѹ��ϵ����a=1-Kd/Ks,��Ks>>Kd���ԣ�a=1.Kd,Ks,Kf�ֱ��ʾ����Ǽܵ����ģ���������������ģ���Լ��������ģ��
k_dynamic=k_darcy/(density_water*g);%��̬������͸ϵ��k_dynamic=k_darcy/(desity_water*g)



for i=1:101%���i��ʾ���������߶� 
    i1=i-1;
  for  j=1:1%���ѭ����ʾ��������ʱ���ѭ��
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
    
fun1=@(s) heaviside(s-i1);%����û����f��tao���൱��f(tao)=1;

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








