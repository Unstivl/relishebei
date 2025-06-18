clear
clc
close all
load('cauResult.mat');

x = 0.5:1:(0.5 + N_p/N_fp);
y = 0.5:1:H/deltaZ+0.5;
[X, Y] = meshgrid(x, y);

for xx = 1:length(x)
    for yy = 1:length(y)
        if (-1)^xx > 0
            zz = (xx-1)*length(y) + yy; 
            z(xx,yy) = zz;
            if zz>length(T_wall_mean)
                T_wall_z(xx,yy) = T_wall_mean(end);
                T_salt_z(xx,yy) = T_salt(end);
            else
                T_wall_z(xx,yy) = T_wall_mean(zz);
                T_salt_z(xx,yy) = T_salt(zz);
            end
        else
            zz = (xx-1)*length(y) + (length(y) - yy + 1);
            z(xx,yy) = zz;
            if zz>length(T_wall_mean)
                T_wall_z(xx,yy) = T_wall_mean(end);
                T_salt_z(xx,yy) = T_salt(end);
            else
                T_wall_z(xx,yy) = T_wall_mean(zz);
                T_salt_z(xx,yy) = T_salt(zz);
            end
        end
    end
end
z=z';
T_wall_z = T_wall_z';
T_salt_z = T_salt_z';


minVal = min(min([T_wall_z(:); T_salt_z(:)]));
maxVal = max(max([T_wall_z(:); T_salt_z(:)]));

figure;
subplot(1, 2, 1);
imagesc(x, y, T_wall_z);
clim([minVal maxVal])
xlabel('角向位置（以panel计）', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('高度（以网格计）', 'FontSize', 14, 'FontWeight', 'bold');
title('HTM模型吸热管壁面温度', 'FontSize', 14, 'FontWeight', 'bold');
axis equal tight;

subplot(1, 2, 2);
imagesc(x, y, T_salt_z);
clim([minVal maxVal])
xlabel('角向位置（以panel计）', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('高度（以网格计）', 'FontSize', 14, 'FontWeight', 'bold');
title('HTM模型吸热管熔融盐温度', 'FontSize', 14, 'FontWeight', 'bold');
axis equal tight;

colormap('autumn');
colorbar;