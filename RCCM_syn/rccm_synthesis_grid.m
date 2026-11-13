
% function [alpha_opt,mu_opt,W_coef,Y_coef,w_poly,y_poly] = rccm_synthesis_grid(plant,controller,W_list, ,Y_list,Y_eval,w_poly,w_poly_fcn,dw_poly_dt_fcn,y_poly,y_poly_fcn,lambda,state_set)
% function       = rccm_syn_grid(plant, rccm ,  )
%%
% adapted from fct: getPeak2peak_gen for paper 'safety-critical control for ensemble systems'

% x is the variables that appear in W and L
c_x1=rccm.c_x1; 
c_u=rccm.c_u; % for mincx c_x*ga_x +c_u*ga_u
options=rccm.options; 
%      options=[1e-4 500 0 10 1];
ccm_eps = rccm.ccm_eps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = plant.n; nu = plant.nu;
nw = plant.nw;
ny_1= size(plant.C,1);
ny_2= size(plant.D,1);
A = plant.A;
B = plant.B;
Bw = plant.Bw;
C = plant.C;
D = plant.D;
nW = length(rccm.W_poly);
nL = length(rccm.L_poly);
W_poly = rccm.W_poly;
L_poly = rccm.L_poly;
range=rccm.al_lim;
Np=rccm.al_np;
f_dyn = plant.dynamics;
if rccm.al_lim(1)==rccm.al_lim(2)
    alpha_temp = range(1);
else
  alpha_temp=linspace(range(1), range(2),Np);
end

alpha_vec=alpha_temp;

%%%%%%%%%%%%%%%(
p2p_vec=1e10.*alpha_vec;

for i= 1: n
    x_range{i} = linspace(rccm.x_lim(i,1),rccm.x_lim(i,2), rccm.x_np(i)); 
end

for i=1:nu
    u_range{i} = linspace(rccm.u_lim(i,1), rccm.u_lim(i,2), rccm.u_np(i));
end

w_range= rccm.w_lim;

d_range = rccm.d_lim;

% Note that LMI1 are affine w.r.t to u; therefore, evaluation at the
% vertices is enough
if nw>1
[T{1}, T{2}, T{3},  T{4}, T{5}, T{6}, T{7}, T{8}]=ndgrid(x_range{1}, x_range{2}, u_range{1}, u_range{2}, u_range{3}, w_range(:,1),w_range(:,2), d_range);
Var_com=zeros(numel(T{1}), length(T));   
else  % nw=1

if nu==3 
    % for 3 inputs
[T{1}, T{2}, T{3},  T{4}, T{5}, T{6}, T{7}]=ndgrid(x_range{1}, x_range{2}, u_range{1}, u_range{2}, u_range{3}, w_range, d_range);
Var_com=zeros(numel(T{1}), length(T));
else
% for 2 inputs
[T{1}, T{2}, T{3},  T{4}, T{5}, T{6}]=ndgrid(x_range{1}, x_range{2}, u_range{1}, u_range{2}, w_range,d_range);
Var_com=zeros(numel(T{1}), length(T));
end
end

for i=1:length(T)
  Var_com(:,i)= reshape(T{i},[],1 );
end



for i=1:length(alpha_vec)
alpha=alpha_vec(i);
 setlmis([]);
 ga=lmivar(1,[1,1]);
 nv=lmivar(1,[1,1]);
 ga_x1 = lmivar(1,[1 1]);  % for x_1
 ga_x2 = lmivar(1, [1 1]);  % for x_2
 ga_u = lmivar(1, [1 1]);

  for k=1:nW+1
  W(k) = lmivar(1,[n,1]);
  end

  for kk=1: nL+1
  L(kk) = lmivar(2,[nu,n]);
  end
    for j=1:size(Var_com,1)
     if nW~=0
        x_1=Var_com(j,1); x_2=Var_com(j,2);  d = Var_com(j,n+nu+nw+1:end); 
       W_coff = [1 subs(W_poly)];
       temp = f_dyn (Var_com(j,n+nu+1:n+nu+nw), Var_com(j,1:n)', Var_com(j,n+1:n+nu)', Var_com(j,n+nu+nw+1:end));
       x_1=temp(1); x_2=temp(2);
       dW_coff = [1 subs(W_poly)];
     else 
         W_coff=1;
         dW_coff=1;
     end
     if nL~=0
         x_1=Var_com(j,1); x_2=Var_com(j,2);  
         d = Var_com(j,n+nu+nw+1:end);
         L_coff = [1 subs(L_poly)];
     else
         L_coff = 1;
     end

      A_p = A(Var_com(j,:));
      B_p = B(Var_com(j,:));
      Bw_p= Bw(Var_com(j,:));
      

     %%%% W>ccm_eps*I
     ni=newlmi;
     for k=1:nW+1
     lmiterm([-ni,1, 1, W(k)], W_coff(k),1);    
     end
     lmiterm([ni, 1,1, 0], ccm_eps);

     %%%% LMI1 for contraction
     ni=newlmi;
     for k=1:nW+1
        
     lmiterm([-ni, 1,1, W(k)], dW_coff(k),1 );
     lmiterm([-ni, 1,2, W(k)], W_coff(k)*A_p, 1 );
     lmiterm([-ni, 2,2, W(k)], (1-alpha)*W_coff(k), 1);
     end
     lmiterm([-ni, 1,3, 0], Bw_p);
     lmiterm([-ni, 3,3 ,nv], 1,1);
     for kk=1:nL+1
         lmiterm([-ni, 1,2, L(kk)], B_p, L_coff(kk));
     end


     %%%% LMI2.1 LMI2.2  for state x_1 x_2 with \ga_x1 \ga_x2
     ni=newlmi;
     for k=1:nW+1
     lmiterm([-ni, 1, 1, W(k)], W_coff(k)*alpha, 1); 
     lmiterm([-ni, 3,1, W(k) ], W_coff(k)*C(1,:),1);
     end
     lmiterm([-ni, 2, 2, ga_x1], 1, 1);
     lmiterm([-ni, 2, 2, nv], -1, 1);
     lmiterm([-ni, 3,3, ga_x1], 1, 1);
    
     ni=newlmi;
     for k=1:nW+1
     lmiterm([-ni, 1, 1, W(k)], W_coff(k)*alpha, 1); 
     lmiterm([-ni, 3,1, W(k) ], W_coff(k)*C(2,:),1);
     end
     lmiterm([-ni, 2, 2, ga_x2], 1, 1);
     lmiterm([-ni, 2, 2, nv], -1, 1);
     lmiterm([-ni, 3,3, ga_x2], 1, 1);

     %%% LMI3 for output with \ga_u
     ni=newlmi;
     for k=1:nW+1
     lmiterm([-ni, 1, 1, W(k)], W_coff(k)*alpha, 1); 
     end
     for kk=1:nL+1
     lmiterm([-ni, 3,1, L(kk) ], L_coff(kk)*D,1);
     end
     lmiterm([-ni, 2, 2, ga_u], 1, 1);
     lmiterm([-ni, 2, 2, nv], -1, 1);
     lmiterm([-ni, 3,3, ga_u], 1, 1);
    end

    lmi=getlmis;
    ndec=decnbr(lmi);
    c=zeros(ndec,1);
    for k=1:ndec
        if defcx(lmi,k,ga_x1)>0; c(k)=c_x1;end
%         if defcx(lmi,k,ga_x2)>0; c(k)=c_x2;end
        if defcx(lmi,k,ga_u)>0; c(k)=c_u; end
    end

    % tic
    [gamma,xopt]=mincx(lmi,c,options);
    % toc
  
disp('-----------------------------------------------------------')
disp(['Optimization progress: ', num2str(i),  '/', num2str(length(alpha_vec))', ' finished!']);



    if ~isempty(xopt)
        p2p_vec(i)=gamma;
        for k=1:nW+1
        W_i{i}(:, n*(k-1)+1:k*(n))=dec2mat(lmi,xopt,W(k));
        end
        for kk=1:nL+1
        L_i{i}(:, n*(kk-1)+1:kk*(n))=dec2mat(lmi,xopt,L(kk));
        end
        nv_i(i)=dec2mat(lmi,xopt,nv);
       ga_x1_i(i)=dec2mat(lmi, xopt, ga_x1);
       ga_x2_i(i)=dec2mat(lmi, xopt, ga_x2);
       ga_u_i(i)=dec2mat(lmi, xopt, ga_u);
   
    else
        disp(['Infeasible at alpha = ', num2str(alpha_vec(i))]);
    end

end
%%
[ga_gain, ind] = min(p2p_vec);
ga_x1_opt = ga_x1_i(ind);
ga_x2_opt = ga_x2_i(ind);
ga_u_opt = ga_u_i(ind);
gain_x1 = sqrt(ga_x1_opt*alpha_vec(ind));
gain_x2 = sqrt(ga_x2_opt*alpha_vec(ind));
gain_u = sqrt( ga_u_opt*alpha_vec(ind));
alpha_opt=alpha_vec(ind);
W_opt = W_i{ind};
L_opt = L_i{ind};
nv_opt =  nv_i(ind);
disp('-----------------------------------------------------------')
disp('-----------------------------------------------------------')
disp(['gain für state x1 : ' num2str(gain_x1), ' at alpha=', num2str(alpha_vec(ind))]);
disp(['gain für state x2  : ' num2str(gain_x2), ' at alpha=', num2str(alpha_vec(ind))]);
disp(['gain für input  : ' num2str(gain_u), ' at alpha=', num2str(alpha_vec(ind))]);
disp('-----------------------------------------------------------')


rccm.ga_x_opt=[ga_x1_opt ga_x2_opt];
rccm.ga_u_opt=ga_u_opt ;
rccm.alpha_opt=alpha_opt;
rccm.W_opt=W_opt;
rccm.L_opt=L_opt; 
rccm.nv_opt=nv_opt;
rccm.alpha_opt= alpha_opt;
 

save(res_name{1}, "plant", "rccm");



