%% Generate gridding points for ioss_synthesis_grid.m
%% Outputs : Var_com (benchmark MHE), Var_com_s (two stages MHE)

for i= 1: length(rccm.x_np)
    x_range{i} = linspace(rccm.x_lim(i,1),rccm.x_lim(i,2), rccm.x_np(i)); 
end



for i=1:length(rccm.u_np)
    u_range{i} = linspace(rccm.u_lim(i,1), rccm.u_lim(i,2), rccm.u_np(i));
end

for i= 1: length(rccm.w_np)
    w_range{i} = linspace(rccm.w_lim(i,1),rccm.w_lim(i,2), rccm.w_np(i)); 
end

for i=1:length(rccm.d_np)
    d_range{i} = linspace(rccm.d_lim(i,1), rccm.d_lim(i,2), rccm.d_np(i));
end


for i= 1: length(rccm.v_np)
    v_range{i} = linspace(rccm.v_lim(i,1),rccm.v_lim(i,2), rccm.v_np(i)); 
end

for i=1:length(rccm.y_np)
    y_range{i} = linspace(rccm.y_lim(i,1), rccm.y_lim(i,2), rccm.y_np(i));
end


for i=1:length(rccm.z_np)
    z_range{i} = linspace(rccm.z_lim(i,1), rccm.z_lim(i,2), rccm.z_np(i));
end



clear T
[T{1}, T{2}, T{3},  T{4}, T{5}, T{6}, T{7}]=ndgrid(x_range{1}, w_range{1}, w_range{2},  d_range{1}, d_range{2}, d_range{3}, u_range{1});
Var_com=zeros(numel(T{1}), length(T));   
for i=1:length(T)
  Var_com(:,i)= reshape(T{i},[],1 );
end

clear T_s
[T_s{1}, T_s{2}, T_s{3},  T_s{4}, T_s{5}, T_s{6}]=ndgrid(z_range{1}, z_range{2}, v_range{1},  v_range{2}, y_range{1}, u_range{1});
Var_com_s=zeros(numel(T_s{1}), length(T_s)); 
for i=1:length(T_s)
  Var_com_s(:,i)= reshape(T_s{i},[],1 );
end

rccm.Var_com=Var_com;
rccm_f.Var_com = Var_com;
rccm_s.Var_com = Var_com_s;
