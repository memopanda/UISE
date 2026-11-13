function [M_d] = MatrixDilation_withZero(M,ind_v,par_r,par_c)
%UNTITLED Summary of this function goes here
% M input matrix
%  e.g. M=[M11 M12 ; M21 M22] \in R^((nx1+nx2)*(ny1+ny2)),  ind_v=[1 0 1]; par_r=[nx1 r1 nx2]
%  par_c=[ny1 r2 ny2];
%  M_d=[M11 O(nx1, r2) M12; O(r1,ny1+r2+ny2); M21 O(nx2,r2) M22];
Rt=zeros(size(M,2), sum(par_c));
Lt=zeros(sum(par_r), size(M,1));
ind1=find(ind_v);
M_par_c=par_c(ind1);
M_par_r=par_r(ind1);
for i=1:sum(ind_v)
    id_Rt_r=sum(M_par_c(1:i))-M_par_c(i)+1:sum(M_par_c(1:i));
    id_Rt_c=sum(par_c(1:ind1(i)))-par_c(ind1(i))+1 :sum(par_c(1:ind1(i)));
    Rt(id_Rt_r,id_Rt_c)=eye(M_par_c(i)); 
    id_Lt_c=sum(M_par_r(1:i))-M_par_r(i)+1:sum(M_par_r(1:i));
    id_Lt_r=sum(par_r(1:ind1(i)))-par_r(ind1(i))+1 :sum(par_r(1:ind1(i)));
    Lt(id_Lt_r,id_Lt_c)=eye(M_par_r(i)); 
end
M_d=Lt*M*Rt;
end