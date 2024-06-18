%��������׼ȷ�ģ�׼ȷ�Լ�Simon����ͼ4
%��������ͼ��ͼ����ֱ����ͼ����ӱ༭��
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




  for  j=1:81 %���j��ʾ����ʱ��
  fprintf('j %f***ѭ������***\r\n',j);      
    
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
xlabel('\tau=t/��k','fontsize',15);
ylabel('$${\hat{u}}=uV_c/k\sigma_0$$','interpreter','latex','fontsize',15,'fontweight','bold')
%('\hat{u}=uV_c/k\sigma_0','interpreter','latex'); 

 