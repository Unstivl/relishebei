clc
clear

d_out = 42.2e-3; % 外径
B = 44.2e-3; % 管圆心距离
Ns = 37; % 柱坐标角向网格数

deltaZ = 13*d_out; % 柱坐标纵向网格长度
deltaPhi = pi/Ns; % 柱坐标角向网格角宽度

sigma = 5.67e-8; % Stefan-Boltzmann constant

Farray = zeros(Ns, Ns+2);

for m = 1:Ns
    for j = 1:Ns+1
        Farray(m,j) = F(m,j,B,d_out,deltaPhi,Ns);
    end
    Farray(m,Ns+2) = F(m,0,B,d_out,deltaPhi,Ns);
end

% Calculate the sum of each row in Farray
rowSums = sum(Farray, 2);
mean(abs(rowSums - 1))



















function val = kroneckerDelta(m, n)
    if m == n
        val = 1;
    else
        val = 0;
    end
end

function val = F(m,j,B,d_out,deltaPhi,Ns) % 第m块微元对第j块微元的角系数
    if j== Ns+1 % 对保温层
        Jd = sqrt( B.^2 + (0.5.*d_out).^2 - B.*d_out.*sin(m.*deltaPhi) );
        alpha_ = acos(d_out./(2*Jd)) - asin(d_out.*cos(m.*deltaPhi))./(2*Jd);
        val = 0.5*(1 - cos(alpha_- pi/2 + m.*deltaPhi));
        if (m+0.5).*deltaPhi <= asin(d_out./B)
            val = 0; % 遮挡
        end
    elseif j== 0 % 对空气
        Jc = sqrt( B.^2 + (0.5.*d_out).^2 - B.*d_out.*sin((m+1).*deltaPhi) );
        alpha_ = acos(d_out./(2*Jc)) + asin(d_out.*cos((m+1).*deltaPhi))./(2*Jc);
        val = 0.5*(1 + cos(pi/2 - alpha_ + (m+1).*deltaPhi));
        if (m+0.5).*deltaPhi >= pi - asin(d_out./B)
            val = 0; % 遮挡
        end
    elseif j < Ns+1 && j > 0
        ac = sqrt( (B - 0.5.*d_out.*(sin((j+1).*deltaPhi) + sin((m+1).*deltaPhi))).^2 + (0.5.*d_out.*(cos((j+1).*deltaPhi) - cos((m+1).*deltaPhi))).^2 );

        bd = sqrt( (B - 0.5.*d_out.*(sin(j.*deltaPhi) + sin(m.*deltaPhi))).^2 + (0.5.*d_out.*(cos(j.*deltaPhi) - cos(m.*deltaPhi))).^2 );

        ad = sqrt( (B - 0.5.*d_out.*(sin((j+1).*deltaPhi) + sin(m.*deltaPhi))).^2 + (0.5.*d_out.*(cos((j+1).*deltaPhi) - cos(m.*deltaPhi))).^2 );

        bc = sqrt( (B - 0.5.*d_out.*(sin(j.*deltaPhi) + sin((m+1).*deltaPhi))).^2 + (0.5.*d_out.*(cos(j.*deltaPhi) - cos((m+1).*deltaPhi))).^2 );

        cd = d_out./2.*deltaPhi;

        val = max( (ad + bc - ac - bd)./(2.*cd) , 0);
    else
        error('j must be in the range [0, Ns+1]');
    end
end