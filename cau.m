clc
clear

sig = 5.67e-8; % Stefan-Boltzmann constant
d_int = 38.9e-3;
d_out = 42.2e-3; % 外径
B = 44.2e-3; % 管圆心距离
H = 7.4; % 吸热器高度
T_salt_0 = 290+273.15; % 初始盐温度
T_sky = 19.5+273.15; % 天空温度
T_amb = 30+273.15; % 环境温度
N_fp = 2; % 吸热器路径数
N_p = 18; % z方向上的panel数（Fig1），是多条路径的网格数的和
N_t = 24; % 吸热器总的网格数（熔融盐流动方向）
epsilon_t = 0.87; % 吸热管发射率
epsilon_sky = 0.895;
epsilon_N_s1 = 0.2;
epsilon_gr = 0.955;
alpha_ = 0.96; % 吸热管吸收率
R_foul = 8.8e-5;
v = 15; % 风速m/s
m_ = 281.6; % 熔融盐流量kg/s，case C


Ns = 37; % 柱坐标角向网格数

deltaZ = 13*d_out; % 柱坐标纵向网格长度
deltaPhi = pi/Ns; % 柱坐标角向网格角宽度

sigma = 5.67e-8; % Stefan-Boltzmann constant

F_ = zeros(Ns+1, Ns+2);

for m = 1:Ns
    for jj = 1:Ns+1
        F_(m,jj) = F(m,jj,B,d_out,deltaPhi,Ns);
    end
    F_(m,Ns+2) = F(m,0,B,d_out,deltaPhi,Ns); % 对环境的角系数
end

% 第Ns+1行是对环境，即0的角系数；这一行第Ns+1列为1减前Ns列，Ns+2列为自身对自身=0
for jj = 1:Ns
    F_(Ns+1,jj) = F0(jj,B,d_out,deltaPhi,Ns);
end

F_(Ns+1,Ns+1) = 1 - sum(F_(Ns+1,1:Ns)); % 对环境的角系数归一性
F_(Ns+1,Ns+2) = 0; % 自己对自己

% 验证角系数归一性
% rowSums = sum(F_, 2);
% mean(abs(rowSums - 1))

%%

T_wall = (290+273.15)*ones(floor(H*N_p/N_fp/deltaZ),Ns); % 壁面温度
T_Ns1 = (290+273.15)*ones(floor(H*N_p/N_fp/deltaZ),1);
q_0 = zeros(floor(H*N_p/N_fp/deltaZ),1); % 出射环境热流
q_j = zeros(floor(H*N_p/N_fp/deltaZ),Ns); % 吸热管壁面辐射热流
q_h = 16e5*ones(floor(H*N_p/N_fp/deltaZ),1); % 镜场热辐射决定的
% for ii = 1:floor(H*N_p/N_fp/deltaZ)
%     z = ii*deltaZ - H*floor(ii*deltaZ/H);
%     q_h(ii) = q_h(ii).*sin(z/H*pi); 
% end % 一种假设的分布：沿着吸热器周向均匀，轴向sin

for ii = 1:floor(H*N_p/N_fp/deltaZ)
    q_h(ii) = q_h(ii)* cos(floor(ii*deltaZ/H)/floor(N_p/N_fp)*pi/2);
end % 一种假设的分布：沿着吸热器轴向均匀，周向cos

T_0 = ((epsilon_sky*T_sky^4 + epsilon_gr*T_amb^4)/(epsilon_gr+epsilon_sky))^0.25; % 环境温度

% for z_num = 1:floor(H*N_p/N_fp/deltaZ)
    % z = (z_num-0.5)*deltaZ; % 取管道上网格的中心
    z_num = 1

    eps0 = epsilon_gr;
    RadA = zeros(Ns+1,Ns+1);
    for m=1:Ns
        RadA(m,1) = (kDelta(m,0)/eps0-(1/eps0 - 1)*F_(m,Ns+2))/sig; % 对应q_0
        for jj=1:Ns
            RadA(m,jj+1) = (kDelta(m,jj)/epsilon_t-(1/epsilon_t - 1)*F_(m,jj))/sig;
        end % 对应q_j
    end % 第一行到第Ns行，代表m=1,...,Ns

    % m=0情况
    RadA(Ns+1,1) = (kDelta(0,0)/eps0-(1/eps0 - 1)*F_(Ns+1,Ns+2))/sig; % 对应q_0
    for jj=1:Ns % 对应q_j
        RadA(Ns+1,jj+1) = (kDelta(0,jj)/epsilon_t-(1/epsilon_t - 1)*F_(Ns+1,jj))/sig;
    end % 最后一行，代表m=Ns+1

    Radb = zeros(Ns+1,1);
    for m=1:Ns
        Radb(m) = (kDelta(m,Ns+1)-F_(m,Ns+1))*T_Ns1(z_num)^4 + (kDelta(m,0)-F_(m,Ns+2))*T_0^4 ;
        for jj=1:Ns
            Radb(m) = (kDelta(m,jj)-F_(m,jj))*T_wall(z_num,jj)^4 + Radb(m);
        end
        Radb(m) = Radb(m) - F_(m,Ns+2)*q_h(z_num)/sig*alpha_;
    end

    Radb(Ns+1) = (kDelta(0,Ns+1)-F_(Ns+1,Ns+1))*T_Ns1(z_num)^4 + (kDelta(0,0)-F_(Ns+1,Ns+2))*T_0^4;
    for jj=1:Ns
        Radb(Ns+1) = (kDelta(0,jj)-F_(Ns+1,jj))*T_wall(z_num,jj)^4 + Radb(Ns+1);
    end
    Radb(Ns+1) = Radb(Ns+1) - F_(Ns+1,Ns+2)*q_h(z_num)/sig*alpha_;

    q = RadA\Radb; % 求解辐射热流
    q_0(z_num) = q(1); % 对应q_0
    q_j(z_num,:) = q(2:Ns+1)'; % 对应q_j
















function val = kDelta(m, n)
    if m == n
        val = 1;
    else
        val = 0;
    end
end



function val = F0(j,B,d_out,deltaPhi,Ns)
    F_0 = F(j,0,B,d_out,deltaPhi,Ns);
    if F_0 <= 0
        val = 0;
    else
        cd = d_out/2 * deltaPhi;
        % ce = sqrt( B.^2 - B.*d_out.*cos((j+1)*deltaPhi-0.5*pi));
        % val = 1/B* (F_0.*cd + ce);
        val = 1/B* (F_0.*cd);
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