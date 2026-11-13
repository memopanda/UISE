%%% eps
 f_w_dt = plant.f_fcn_w_dt;
 if plant.nw==1
F_w = @(z)(-(1e15*f_w_dt(z(1),[z(2); z(3)])-1e15*z(1))^2);  % x_d x_c x_T
 else  % two slow states
F_w = @(z)(-norm(1e15*f_w_dt([z(1);  z(2)],[z(3) ; z(4)])-1e15*[z(1);z(2)]));  % x_d x_c x_T
 end
lb = [ rccm.w_lim(:,1)' rccm.x_lim(:,1)' ];
ub = [ rccm.w_lim(:,end)' rccm.x_lim(:,2)' ];
 a = 0.01:0.01: 0.99;

% options = optimoptions(@fmincon, 'Algorithm', 'sqp');
options = optimoptions('fmincon');
% options.Algorithm = "sqp";
options.StepTolerance = 1e-15;
options.OptimalityTolerance = 1e-15;
options.ConstraintTolerance = 1e-10;
options.Display = 'none';

% val_opt=ones(1,length(a));

for i=1:length(a)
    x0 = a(i)*lb+(1-a(i))*ub;
[x, fval,exitflag, output] = fmincon(F_w, x0, [], [], [],[], lb, ub, [], options);
   x_opt(:, i) = x;

   val_opt(i) = fval;
end
 [val, id]=min(val_opt);
 xx_opt= x_opt(:,id);
%  eps_test= abs(f_w_dt(0.1998, [0.0068 33.3892]) -0.1998);
 eps= norm(f_w_dt(xx_opt(1:plant.nw), xx_opt(plant.nw+1:end)) -xx_opt(1:plant.nw));

%%%% eps_b
 f_w_DT = plant.f_fcn_w_DT;
F_w = @(z)(-norm(1e15*f_w_DT([z(1);  z(2)],[z(3) ; z(4)])-1e15*[z(1);z(2)]));
for i=1:length(a)
    x0 = a(i)*lb+(1-a(i))*ub;
[x, fval,exitflag, output] = fmincon(F_w, x0, [], [], [],[], lb, ub, [], options);
   x_opt(:, i) = x;

   val_opt(i) = fval;
end
 [val, id]=min(val_opt);
 xx_opt= x_opt(:,id);
%  eps_test= abs(f_w_dt(0.1998, [0.0068 33.3892]) -0.1998);
 eps_b= norm(f_w_DT(xx_opt(1:plant.nw), xx_opt(plant.nw+1:end)) -xx_opt(1:plant.nw));

%%% r= Lg_DT* mH_DT

%   fg= plant.dg_fcn_w_DT;
%   dg_w = @(w)(-norm(fg(w)));
% lb = [ rccm.w_lim(:,1)' ];
% ub = [ rccm.w_lim(:, end)'];
% a = 0.01:0.01: 0.99;
% for i=1:length(a)   
%     x0 = a(i)*lb+(1-a(i))*ub;
% [x, fval,exitflag, output] = fmincon(dg_w, x0, [], [], [],[], lb, ub, [], options);
% 
%    w_opt(:, i) = x;
%    val_opt(i) = fval;
% end
%  [val, id]=min(val_opt);
%  xx_opt= w_opt(:,id);
%  Lg_DT= norm(fg(xx_opt'));
%  
%  H_DT= plant.H_fcn_w_DT;
%  H_fcn_DT = @(x)(-norm(H_DT(x)));
%  lb = [rccm.x_lim(:,1)' ];
%  ub = [ rccm.x_lim(:,2)' ];
%  a = 0.01:0.01: 0.99;
%  for i=1:length(a)
%     x0 = a(i)*lb+(1-a(i))*ub;
% [x, fval,exitflag, output] = fmincon(H_fcn_DT, x0, [], [], [],[], lb, ub, [], options);
%    x_opt2(:, i) = x;
%    val_opt(i) = fval;
%  end
%   [val, id]=min(val_opt);
%  xx_opt= x_opt2(:,id);
%  mH_DT= norm(H_DT(xx_opt));
%  r= Lg_DT*mH_DT;

% tube.r=r;
tube.eps=eps;
% tube.nom_w= max([norm(rccm.w_lim(:,1)), norm(rccm.w_lim(:,end))]);
tube.eps_b = eps_b;

save(res_name{1}, "tube", "-append");


