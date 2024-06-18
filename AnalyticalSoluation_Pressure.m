
%�����Ӧ��������"AN ANALYTICAL SOLUTION FOR THE TRANSIENT RESPONSE OF SATURATED
%POROUS ELASTIC SOLIDS"��ͼ2(b)��

%ͼ����ֱ����ͼ�ϱ�ע��

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




  for  j=1:81%���j��ʾ����ʱ��
  
    fprintf('j %f***ѭ������***\r\n',j);  
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
digits(6);%���ڶ���vpa�ľ��ȣ�4��ʾС�������λ

fun=@(s)(k_darcy/Vc).*exp(-b*s/(2*a)).*besseli(1,b*sqrt(s.^2-a*i.^2)/(2*a)).*heaviside(s-i*sqrt(a))./sqrt(s.^2-a*i.^2);

W_e(j)=-1*(temp1*temp2+temp1*temp3*vpa(int(fun,s,0.0,tao)));%-1�ǿ�������f(t)=-1;
    

u_e(j)=-(-1)*(k_darcy/Vc)*heaviside(tao-i)-beta*W_e(j);

p(j)=Vc*(beta*u_e(j)+kappa*W_e(j))/k_darcy;

 if(p(j)>=1)%%�������������Ĵ�����Ϊ��������Щ������
     p(j)=1.0;
 end


 end

p

figure (2)
%plot(0.0069:2.0078:200,p,'--r','linewidth', 1.5)
plot(0:2.5:200,p,'--r','linewidth', 2.5)
xlabel('\tau=t/��k','fontsize',15);
 ylabel('$${\hat{p}}=p/\sigma_0$$','interpreter','latex','fontsize',15,'fontweight','bold')
gca=legend({'Analytical Solution'},'FontSize',14,'Location','SouthEast');
po=get(gca,'Position');%��ȡlegend��λ��
set(gca,'Position',[po(1)-0.05,po(2)+0.2,po(3),po(4)]);

hold on
