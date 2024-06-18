clear all
clc

tic;
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
tao=200;%tao不能大于400，超过500之后就不是很吻合，解析解中的无限边界，而我们算的是固定边界，所以会有区别，见宋佳博士论文107页

totlstep=tao*k_darcy*density_solid/(10^-4);%这里的时间t要与Simon文章中的时间tao进行换算，比如，这里的t=0.0598对应的是Simon文章中的tao=40；


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

% Four corner points 四个角点坐标
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
botNodes   = unique(botEdge); %unique函数——取集合的单值元素 用法：b=unique(a)，取集合a的不重复元素构成的向量。
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
% 下面的程序是集成总刚 "K" 和总质量矩阵 "M" %
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
diag0=sum(M,2);%这里将M和S进行对角化的目的是为了求逆方便，如果matlab允许，不对角化进行直接求逆是符合计算公式的
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
%                    定义初值                   %
%乡6%%%%%96%%%%%%%%%%%%%%%%%乡6%%%%%%%%%%%%%%%%% 
fprintf('定义初值');
Fp=zeros(totln,1); 
Fpnext=zeros(totln,1);
fu=zeros(2*totln,1); 
funext=zeros(2*totln,1); 
us_before=zeros(2*totln,1);
us=zeros(2*totln,1); 
usnext=zeros(2*totln,1);
vs=zeros(2*totln,1); 
vsnext=zeros(2*totln,1);
vvs=zeros(2*totln,1);%初始的加速度
pressure=0.0*ones(totln,1); 
pressurenext=zeros(totln,1);
%uw=zeros(2*totln,1); 
%uwnext=zeros(2*tofln,1); 
%ur=textread('dizhenweiyi?kobe.txt’');

 theta1=1;
 theta2=1;
%%乡5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                节点载荷施加                    %
%%%多6%%%%%%%%%%%%%%%%%%%%%%%多5%%%%%%%%%%%%%%%%%% 
   for i=1:totln                   
        if(i==201||i==202)  %在节点11和22上      
           fu(2*i)=-1.0*(zonewidth/XmeshNum/2);   %对于一维的问题，认为这里1的载荷全部加在一个节点上，不能除以2     
        end                                                %加载的力要与Simon文章中的f(t)要进行换算，加载的力换算了之后，最后计算的位移就不需要换算了
   end 


%%乡5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                  逐步迭代计算                   %
%%%多6%%%%%%%%%%%%%%%%%%%%%%%多5%%%%%%%%%%%%%%%%%% 
for tt=1:totlstep
  time=tt*dt; 
 fprintf('time %f***开始逐步迭代计算***\r\n',time);
 
   for i=1:totln
      if((i==201)||(i==202))
       pressurenext(i)=0.0;
       pressure(i)=0.0;
      end
    end 

   % us_before=us-vs*dt+0.5*vvs*dt^2;
   
 
%%%%%%%%%%%%%%%%%%%%%%位移计算%%%%%%%%%%%%%%%%%%%% 

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
%%%%%%%%%%透射边界%%%%%%%%%%位移边界条件

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

   %%%%%%%%%%%%%%%%%%%%%%孔压计算%%%%%%%%%%%%%%%%%%%%
   temp2=S+theta1*dt*J;
   temp3=Fp-J*pressure-Q'*vs;
   B_n=inv(temp2)*temp3; 
   
   
   
  pressurenext=pressure+B_n*theta1*dt;

  % pressurenext=pressure+dt*inv(S)*(-J*pressure-Q'*vsnext);
 
         for i=1:2*totln
            us(i)=usnext(i);
            vs(i)=vsnext(i); 
         end
%%%%%%%%%%孔压边界条件%%%%%%%%%% 
 for i=1:totln
   if((i==201)||(i==202))
       pressurenext(i)=0.0;
       pressure(i)=0.0;
   end
 end 
 
 

  for i=1:totln
      pressure(i)=pressurenext(i); 
  end
  

 %%%%%乡6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %                  输出信息                    %
 
  if(tt==1||mod(tt,37)==0) %mod是取余运算
  %  fprintf(fid1,'%f %f\r\n',time/(0.004883*density_solid),us(402)*Vc/0.004883); 
   fprintf(fid1,'%f %f\r\n',time/(0.004883*density_solid),usnext(402)*Vc/k_darcy); %获得顶点在0到tao的位移变化
  fprintf(fid2,'%f %f\r\n',time/(0.004883*density_solid),pressurenext(193)); 
  end
  temp=usnext(2:4:402);%这里是另一个算例，想用1到100个节点在tao时刻的位移与解析解的结果，发现吻合不是太好。
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
xlabel('\tau=t/ρk','fontsize',15);
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
 po=get(gca,'Position');%获取legend的位置
 set(gca,'Position',[po(1)-0.05,po(2)+0.3,po(3),po(4)]);
% plot(time1,scpp3,'g','linewidth', 1.5)
% hold on
% plot(time1,scpp4,'g','linewidth', 1.5)
% hold on
 xlabel('\tau=t/ρk','fontsize',15);
 ylabel('$${\hat{p}}=p/\sigma_0$$','interpreter','latex','fontsize',15,'fontweight','bold')
 hold on



tic;
toc;