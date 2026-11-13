function [x_est, err]=proj_zx(z, ym, T_mx, x_set)
% project z onto the set \X
% output: x in \X
% Only work for  T_mx*([ym(1)+\eps; z]), T_mx=(1 0 0 ;0 1 0; -be 0 1) or more general: dim of x -dim of z =1
be = -T_mx(3,1);
x_1_range= [x_set(1,1)-ym(1), x_set(1,2)- ym(1) ];
x_3_range =(1/be)*[z(2)-x_set(3,2)-be*ym(1),  z(2)-x_set(3,1)-be*ym(1)];

if (x_1_range(2)<x_3_range(1)) || (x_3_range(2)<x_1_range(1))
    disp('No feasible solution!');
    x_est = [];
else
   x_range = [x_1_range, x_3_range];
   [~, ind] = min(x_range);
   x_range(ind) = [];
   [~,ind2] = max(x_range);
   x_range(ind2) =[];
   ab = x_range;
   if ab(1)*ab(2)<0
       err=0;
       x_est = T_mx*[ym(1);z];
   else
       err= sign(ab(1))*min(abs(ab));
        x_est = T_mx*[ym(1)+sign(ab(1))*min(abs(ab));z];
   end

end

end