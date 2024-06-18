function element=make_elem(node_pattern,num_u,num_v,inc_u,inc_v)

% function element=make_elem(node_pattern,num_u,num_v,inc_u,inc_v)
%
% creates a connectivity list

if ( nargin < 5 )
   disp(['Not enough parameters specified for make_elem function'])
end

inc=[zeros(1,size(node_pattern,2))];%size(node_patern,2)=4
e=1;
element=zeros(num_u*num_v,size(node_pattern,2));%rank=(num_u*num_v)x(4)

for row=1:num_v%row index follow y axis (row)
   for col=1:num_u%row index follow x axis (column)
      element(e,:)=node_pattern+inc;
      inc=inc+inc_u;%inc_u=1=[1 1 1 1]
      e=e+1;%move to next element on the same row(order number of 4 nodes increase by 1)
   end
   inc=row*inc_v;%move from below row to above row (increasing step = 1)
end
