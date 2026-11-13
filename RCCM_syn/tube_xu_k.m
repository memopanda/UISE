function [e_xu, e_w_k]= tube_xu_k(k, es_0, ef_0,  tube, plant, rccm)
bn=plant.DT/plant.dt;
% n0=t-floor(t/bn)*bn;
% e_w=@(k,n0, e_0)(k*tube.eps+e_0+((1+tube.r)^floor(k+n0/bn)-1)*tube.eps_b);
e_w=@(k, e_0)(k*tube.eps+e_0+floor(k/bn)*tube.eps_b);
n=plant.n;
alpha = rccm.alpha_opt;
ga_x=rccm.ga_x_opt;
ga_u = rccm.ga_u_opt;
nv  = rccm.nv_opt;

e_xu=alpha*(1-alpha)^k*ef_0.*ones(1, n+length(ga_u))*diag([ga_x, ga_u]);


if k>0
for j=1:k
  e_xu=e_xu+nv*alpha*(1-alpha)^(k-j)*e_w(j-1, es_0)^2*ones(1, n+length(ga_u))*diag([ga_x, ga_u]);
end
end

e_xu=sqrt(e_xu)';
e_w_k=e_w(k, es_0);
end
