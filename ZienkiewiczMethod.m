clear all
clc

tic;
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

k_dynamic=0.004883;

lamda=833.3;
alpha=1;
Mu=1250;
Vc=sqrt((lamda+2*Mu+alpha^2*Qb)/density_solid);

Xgrid=2;        
Ygrid=101; 
XmeshNum=1;         
YmeshNum=100; 
zonewidth=2.0;      
zoneheight=100.0; 
dt=10^-4; 
tao=200;%tao���ܴ���400������500֮��Ͳ��Ǻ��Ǻϣ��������е����ޱ߽磬����������ǹ̶��߽磬���Ի������𣬼��μѲ�ʿ����107ҳ

totlstep=tao*k_darcy*density_solid/(10^-4);%�����ʱ��tҪ��Simon�����е�ʱ��tao���л��㣬���磬�����t=0.0598��Ӧ����Simon�����е�tao=40��


% +++++++++++++++++++++++++++++
%            MESHING
% +++++++++++++++++++++++++++++
% ---------------------------------------
disp('MESH GENERATION')
% Number of nodes along two directions
L = 2 ;
D = 100 ;
nnx = 2 ;
nny = 101 ;

% Four corner points �ĸ��ǵ�����
pt1 = [0 0] ;
pt2 = [L 0] ;
pt3 = [L D] ;
pt4 = [0 D] ;

% Uniform meshing with Q4 elements
% Data structures for nodes and elements
% Node = [x1 y1
%         x2 y2
%         ...
%         xn yn]
% element = [1 3   5 7
%            4 20 35 78
%            ...
%           ]
elemType = 'Q4' ; 
[node,element] = meshRectangularRegion(pt1, pt2, pt3, pt4, nnx,nny,elemType);

% compute number of nodes, of elements
numnode = size(node,1);
numelem = size(element,1);


% define essential boundaries
uln = nnx*(nny-1)+1;       % upper left node number61.5
urn = nnx*nny;             % upper right node number
lrn = nnx;                 % lower right node number
lln = 1;                   % lower left node number
cln = nnx*(nny-1)/2+1;     % node number at (0,0)

topEdge  = [ uln:1:(urn-1); (uln+1):1:urn ]';
botEdge  = [ lln:1:(lrn-1); (lln+1):1:lrn ]';

% GET NODES ON DIRICHLET BOUNDARY AND ESSENTIAL BOUNDARY
botNodes   = unique(botEdge); %unique��������ȡ���ϵĵ�ֵԪ�� �÷���b=unique(a)��ȡ����a�Ĳ��ظ�Ԫ�ع��ɵ�������
topNodes   = unique(topEdge);
leftNodes=find(node(1:size(node,1),1)<=0.01)
rightNodes=find(abs(node(1:size(node,1),1)-2)<=0.01);


rightEdge=[ rightNodes(1:size(rightNodes,1)-1)'; rightNodes(2:size(rightNodes,1))']';
leftEdge=[ leftNodes(1:size(rightNodes,1)-1)'; leftNodes(2:size(rightNodes,1))']';
dispNodes = [botNodes];
tracNodes = topNodes;

totln=numnode; 
% Plot mesh and enriched nodes to check
%figure (1)
%plot_mesh(node,element,elemType,'k-');


% %%%%%%%%%%%%%%%%%%%-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   MESH GENERATION FINISH'])
% Number of nodes along two directions%%%%%%%%%%%%%%

%error('   MESH GENERATION FINISH')


K=zeros(numnode*2,numnode*2);
M=zeros(numnode*2,numnode*2);
Q=zeros(numnode*2,numnode);
J=zeros(numnode,numnode);
S=zeros(numnode,numnode);
invmass0=zeros(numnode*2,numnode*2);
invdamp0=zeros(numnode,numnode);
[number_element,n1]=size(element);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����ĳ����Ǽ����ܸ� "K" ������������ "M" %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:number_element
    x1=node(element(i,1),1);
    y1=node(element(i,1),2);
    x2=node(element(i,2),1);
    y2=node(element(i,2),2);
    x3=node(element(i,3),1);
    y3=node(element(i,3),2);
    x4=node(element(i,4),1);
    y4=node(element(i,4),2);   
  k=BilinearQuadElementStiffness(E,NU,h,x1,y1,x2,y2,x3,y3,x4,y4,1);
  m=BilinearQuadElementMass(h,density_solid,x1,y1,x2,y2,x3,y3,x4,y4);
  q=BilinearQuadElementCoplingMatrixNew(h,x1,y1,x2,y2,x3,y3,x4,y4);
  jmatrix=BilinearQuadElementConductiveMatrix(h,k_dynamic,x1,y1,x2,y2,x3,y3,x4,y4);
  s=BilinearQuadElementCompressMatrix(h,Qb,x1,y1,x2,y2,x3,y3,x4,y4);
  
  K=BilinearQuadAssemble_K_M(K,k,element(i,1),element(i,2),element(i,3),element(i,4));
  M=BilinearQuadAssemble_K_M(M,m,element(i,1),element(i,2),element(i,3),element(i,4));
  Q=BilinearQuadAssemble_Q(Q,q,element(i,1),element(i,2),element(i,3),element(i,4));
  J=BilinearQuadAssemble_J_S(J,jmatrix,element(i,1),element(i,2),element(i,3),element(i,4));
  S=BilinearQuadAssemble_J_S(S,s,element(i,1),element(i,2),element(i,3),element(i,4));
  
end
Q=-Q;
diag0=sum(M,2);%���ｫM��S���жԽǻ���Ŀ����Ϊ�����淽�㣬���matlab�������Խǻ�����ֱ�������Ƿ��ϼ��㹫ʽ��
diag1=sum(S,2);
for i=1:2*totln
    for j=1:2*totln
     if(i==j)
    invmass0(i,j)= diag0(i) ;
     else
    invmass0(i,j) =0; 
     end
    end
end

invmass=inv(invmass0);

for i=1:totln
    for j=1:totln
     if(i==j)
    invdamp0(i,j)= diag1(i) ;
     else
    invdamp0(i,j)=0; 
     end
    end
end
invdamp=inv(invdamp0);

fid1=fopen('displacement.txt','w'); 
fid2=fopen('pressure.txt','w');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                    �����ֵ                   %
%��6%%%%%96%%%%%%%%%%%%%%%%%��6%%%%%%%%%%%%%%%%% 
fprintf('�����ֵ');
Fp=zeros(totln,1); 
Fpnext=zeros(totln,1);
fu=zeros(2*totln,1); 
funext=zeros(2*totln,1); 
us_before=zeros(2*totln,1);
us=zeros(2*totln,1); 
usnext=zeros(2*totln,1);
vs=zeros(2*totln,1); 
vsnext=zeros(2*totln,1);
vvs=zeros(2*totln,1);%��ʼ�ļ��ٶ�
pressure=0.0*ones(totln,1); 
pressurenext=zeros(totln,1);
%uw=zeros(2*totln,1); 
%uwnext=zeros(2*tofln,1); 
%ur=textread('dizhenweiyi?kobe.txt��');

 theta1=1;
 theta2=1;
%%��5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                �ڵ��غ�ʩ��                    %
%%%��6%%%%%%%%%%%%%%%%%%%%%%%��5%%%%%%%%%%%%%%%%%% 
   for i=1:totln                   
        if(i==201||i==202)  %�ڽڵ�11��22��      
           fu(2*i)=-1.0*(zonewidth/XmeshNum/2);   %����һά�����⣬��Ϊ����1���غ�ȫ������һ���ڵ��ϣ����ܳ���2     
        end                                                %���ص���Ҫ��Simon�����е�f(t)Ҫ���л��㣬���ص���������֮���������λ�ƾͲ���Ҫ������
   end 


%%��5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                  �𲽵�������                   %
%%%��6%%%%%%%%%%%%%%%%%%%%%%%��5%%%%%%%%%%%%%%%%%% 
for tt=1:totlstep
  time=tt*dt; 
 fprintf('time %f***��ʼ�𲽵�������***\r\n',time);
 
   for i=1:totln
      if((i==201)||(i==202))
       pressurenext(i)=0.0;
       pressure(i)=0.0;
      end
    end 

   % us_before=us-vs*dt+0.5*vvs*dt^2;
   
 
%%%%%%%%%%%%%%%%%%%%%%λ�Ƽ���%%%%%%%%%%%%%%%%%%%% 

K0=K+Q*inv(S)*Q';

temp1=M+dt^2*theta2*K0/2;

a_n=inv(temp1)*(fu-K0*(us+vs*dt*theta1)+Q*inv(S)*Q'*us+Q*pressure);

vsnext=vs+a_n*theta1*dt;


  for i=1:totln 
           if (i==1||i==2)      
             vsnext(2*i-1)=0;
             vsnext(2*i)=0;  
           else  
             vsnext(2*i-1)=0;
           end  
   end 

usnext=us+vs*theta1*dt+(a_n*theta2*dt^2)/2;
%%%%%%%%%%͸��߽�%%%%%%%%%%λ�Ʊ߽�����

   for i=1:totln 
           if (i==1||i==2)    
              usnext(2*i-1)=0; 
              usnext(2*i)=0;                   
           else  
             usnext(2*i-1)=0; 
           end  
   end 
  
  kama=16/49;
  beta=9/14;
  vsnext=kama/(beta*dt)*(usnext-us)+(1-kama/beta)*vs;
  
   for i=1:totln 
           if (i==1||i==2)    
               vsnext(2*i-1)=0; 
              vsnext(2*i)=0; 
           else         
               vsnext(2*i-1)=0; 
         
           end  
   end 

   %%%%%%%%%%%%%%%%%%%%%%��ѹ����%%%%%%%%%%%%%%%%%%%%
   temp2=S+theta1*dt*J;
   temp3=Fp-J*pressure-Q'*vs;
   B_n=inv(temp2)*temp3; 
   
   
   
  pressurenext=pressure+B_n*theta1*dt;

  % pressurenext=pressure+dt*inv(S)*(-J*pressure-Q'*vsnext);
 
         for i=1:2*totln
            us(i)=usnext(i);
            vs(i)=vsnext(i); 
         end
%%%%%%%%%%��ѹ�߽�����%%%%%%%%%% 
 for i=1:totln
   if((i==201)||(i==202))
       pressurenext(i)=0.0;
       pressure(i)=0.0;
   end
 end 
 
 

  for i=1:totln
      pressure(i)=pressurenext(i); 
  end
  

 %%%%%��6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %                  �����Ϣ                    %
 
  if(tt==1||mod(tt,37)==0) %mod��ȡ������
  %  fprintf(fid1,'%f %f\r\n',time/(0.004883*density_solid),us(402)*Vc/0.004883); 
   fprintf(fid1,'%f %f\r\n',time/(0.004883*density_solid),usnext(402)*Vc/k_darcy); %��ö�����0��tao��λ�Ʊ仯
  fprintf(fid2,'%f %f\r\n',time/(0.004883*density_solid),pressurenext(193)); 
  end
  temp=usnext(2:4:402);%��������һ������������1��100���ڵ���taoʱ�̵�λ���������Ľ���������Ǻϲ���̫�á�
  u_tao=temp(end:-1:1);
end

size(u_tao)

%scvs=textread('velocity.txt');

scus=textread('displacement.txt');  
scpp=textread('pressure.txt'); 

time=scus(:,1); 
scus2=scus(:,2); 
scus2'
% scus3=scus(:,3); 
% scus4=scus(:,4); 
 figure (1)
 plot(time,scus2,'-<b','linewidth', 1.5)
% hold on
% plot(time,scus3,'g','linewidth', 1.5)
% hold on
% plot(time,scus4,'g','linewidth', 1.5)
hold on
xlabel('\tau=t/��k','fontsize',15);
ylabel('$${\hat{u}}=uV_c/k\sigma_0$$','interpreter','latex','fontsize',15,'fontweight','bold')

%figure (1)
%plot(0:1:100,u_tao,'b','linewidth', 1.5)
 
 figure (2)
  time1=scpp(:,1); 
 scpp2=scpp(:,2); 
 
% scpp3=scpp(:,3); 
% scpp4=scpp(:,4); 
% 
 plot(time1,scpp2,'-<b','linewidth', 1.5)
 hold on
 gca=legend({'Zienkiewicz^\primes Method'},'location','SouthEast','FontSize',14);
 po=get(gca,'Position');%��ȡlegend��λ��
 set(gca,'Position',[po(1)-0.05,po(2)+0.3,po(3),po(4)]);
% plot(time1,scpp3,'g','linewidth', 1.5)
% hold on
% plot(time1,scpp4,'g','linewidth', 1.5)
% hold on
 xlabel('\tau=t/��k','fontsize',15);
 ylabel('$${\hat{p}}=p/\sigma_0$$','interpreter','latex','fontsize',15,'fontweight','bold')
 hold on



tic;
toc;