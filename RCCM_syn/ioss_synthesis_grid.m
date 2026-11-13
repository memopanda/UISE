function [rccm] = ioss_synthesis_grid(rccm)
% 
% update on 03.06.2025: extend C, D to C_p, D_p
%%

options=rccm.options; 
%      options=[1e-4 500 0 10 1];
P_eps = rccm.P_eps;
cvy_eps = rccm.cvy_eps;
Var_com = rccm.Var_com;
n = rccm.n;
A = rccm.A;
B = rccm.B;
C = rccm.C;
D = rccm.D;
alpha_vec=linspace(rccm.al_lim(1), rccm.al_lim(2), rccm.al_np);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cx_l_vec=0.1*P_eps*ones(1,length(alpha_vec));
cx_u_vec=ones(1, length(alpha_vec));

feas_sign = 1;

for i=1:length(alpha_vec)
alpha=alpha_vec(i);
 
setlmis([]);

 cx_l=lmivar(1,[1,1]);
 cx_u = lmivar(1,[1,1]);
 cy=lmivar(1,[1,1]);
 cv=lmivar(1,[1,1]);

  P = lmivar(1,[n,1]);

  %%%% cx_u*I > P>cx_l*I, cx_l>P_eps, cvy_eps > cy cv cx_u >1e-3*P_eps, 
     ni=newlmi;
     
     lmiterm([-ni,1, 1, P], 1,1);    
     lmiterm([ni, 1,1, cx_l], 1,1);
     
     ni=newlmi;
     lmiterm([-ni,1, 1, cx_l], 1,1);    
     lmiterm([ni, 1,1, 0], P_eps);
      ni=newlmi;
     lmiterm([-ni,1, 1, cv], 1,1);    
     lmiterm([ni, 1,1, 0], 1e-3*P_eps);
      ni=newlmi;
     lmiterm([-ni,1, 1, cy], 1,1);    
     lmiterm([ni, 1,1, 0], 1e-3*P_eps);

     ni=newlmi;
     lmiterm([ni,1, 1, P], 1,1);    
     lmiterm([-ni, 1,1, cx_u], 1,1);

      ni=newlmi;
     lmiterm([ni,1, 1, cx_u], 1,1);    
     lmiterm([-ni, 1,1, 0], cvy_eps);
     ni=newlmi;
     lmiterm([ni,1, 1, cv], 1,1);    
     lmiterm([-ni, 1,1, 0], cvy_eps);
      ni=newlmi;
     lmiterm([ni,1, 1, cy], 1,1);    
     lmiterm([-ni, 1,1, 0], cvy_eps);

    for j=1:size(Var_com,1)

      A_p = A(Var_com(j,:));
      B_p = B(Var_com(j,:));
      if isa(C,'function_handle') && isa(D,'function_handle')
      C_p = C(Var_com(j,:));
      D_p = D(Var_com(j,:));
      else
          C_p = C;
          D_p = D;
      end
     %%%% LMI1 for contraction
     ni=newlmi;
     
        
     lmiterm([-ni, 1,1, P], alpha,1);
     lmiterm([-ni, 1,1, cy], C_p',C_p);
     lmiterm([-ni, 1,2, cy], 1, C_p'*D_p );
     lmiterm([-ni, 2,2, cv], 1, 1);
     lmiterm([-ni, 2,2, cy], 1, D_p'*D_p);

     lmiterm([-ni, 1,3, P], -A_p',1);
     lmiterm([-ni, 2,3, P], -B_p', 1);
     lmiterm([-ni, 3,3, P], 1,1);
  
    end

    lmi=getlmis;
    ndec=decnbr(lmi);
    c=zeros(ndec,1);
    for k=1:ndec
        if defcx(lmi,k,cx_u)>0; c(k)=1;end
%         if defcx(lmi,k,ga_x2)>0; c(k)=c_x2;end
        if defcx(lmi,k,cx_l)>0; c(k)=-1; end
    end

    % tic
    [~,xopt]=mincx(lmi,c,options);
    % toc
  
disp('-----------------------------------------------------------')
disp(['Optimization progress: ', num2str(i),  '/', num2str(length(alpha_vec))', ' finished!']);

    if ~isempty(xopt)
        P_vec{i}=dec2mat(lmi,xopt,P);
        cx_l_vec(i)=dec2mat(lmi,xopt,cx_l);
        cx_u_vec(i)=dec2mat(lmi, xopt, cx_u);
         cy_vec(i)=dec2mat(lmi, xopt, cy);
         cv_vec(i)=dec2mat(lmi, xopt, cv);
         feas_sign = 0;
   
    else
        disp(['Infeasible at alpha = ', num2str(alpha_vec(i))]);
    end

end
%%
if feas_sign
disp('-----------------------------------------------------------')
disp('-----------------------------------------------------------')
disp('No any feasible solution!');
disp('-----------------------------------------------------------')
else
 cx_vec = cx_u_vec./cx_l_vec;
N_vec = -log(cx_vec*4)./log(alpha_vec);
[N_min, ind] = min(N_vec);
rccm.cx_2_opt = cx_u_vec(ind);
rccm.cx_1_opt = cx_l_vec(ind);
rccm.cy_opt = cy_vec(ind)/alpha_vec(ind);
rccm.cv_opt = cv_vec(ind)/alpha_vec(ind);
rccm.P_opt = P_vec{ind};
rccm.N_min=N_min;
rccm.cx_opt = cx_vec(ind);
rccm.alpha_opt = alpha_vec(ind); 
disp('-----------------------------------------------------------')
disp('-----------------------------------------------------------')
disp(['Minimum Horizon : ' num2str(N_min), ' at alpha=', num2str(alpha_vec(ind))]);
disp('-----------------------------------------------------------')



end

end