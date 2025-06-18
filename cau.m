clc
clear

TOL = 1e-3;
nnMax = 100;
sig = 5.67e-8; % Stefan-Boltzmann constant
d_int = 38.9e-3;
d_out = 42.2e-3; % 外径
B = 44.2e-3; % 管圆心距离
H = 7.4; % 吸热器高度
D = 6; % 吸热器直径
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
k_t = 16;
v = 15; % 风速m/s
m_ = 281.6; % 熔融盐流量kg/s，case C

% 以下计算对流换热
g = 9.81; % 重力加速度
nu = refpropm('V','T',T_amb,'P',1e5,'AIR.MIX')/refpropm('D','T',T_amb,'P',1e5,'AIR.MIX'); % 运动粘度
Pr = refpropm('C','T',T_amb,'P',1e5,'AIR.MIX')*refpropm('V','T',T_amb,'P',1e5,'AIR.MIX')/refpropm('L','T',T_amb,'P',1e5,'AIR.MIX'); % Prandtl number
Gr_H = g*(1/T_amb)*(T_salt_0 - T_amb)*H^3/(nu^2); % Grashof number
Nu_H = 0.11*(Gr_H*Pr)^(1/3); % Nusselt number
h_nc = Nu_H*refpropm('L','T',T_amb,'P',1e5,'AIR.MIX')/H; % 自然对流换热系数
Re_D = v*D/nu; % Reynolds number
Nu_D = 0.18*Re_D^0.63;
h_fc = Nu_D*refpropm('L','T',T_amb,'P',1e5,'AIR.MIX')/D; % 强制对流换热系数
h = (h_nc^3.2 + h_fc^3.2)^(1/3.2);



Ns = 37; % 柱坐标角向网格数

% deltaZ = 13*d_out; % 柱坐标纵向网格长度
deltaZ = H/14; % 柱坐标纵向网格长度
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

T_wall = (10+290+273.15)*ones(floor(H*N_p/N_fp/deltaZ),Ns); % 壁面温度
T_wall_mean = zeros(floor(H*N_p/N_fp/deltaZ),1); % 壁面平均温度
T_Ns1 = (10+290+273.15)*ones(floor(H*N_p/N_fp/deltaZ),1); % 保温层温度
T_salt = (290+273.15)*ones(floor(H*N_p/N_fp/deltaZ),1); % 盐温度
q_0 = zeros(floor(H*N_p/N_fp/deltaZ),1); % 出射环境热流
q_c_l = zeros(floor(H*N_p/N_fp/deltaZ),Ns); % 吸热管壁面对流热流
q_r_l = zeros(floor(H*N_p/N_fp/deltaZ),Ns); % 吸热管壁面辐射热流
q_j = zeros(floor(H*N_p/N_fp/deltaZ),Ns); % 吸热管壁面辐射热流
q_h = 1.6e6*ones(floor(H*N_p/N_fp/deltaZ),1); % 镜场热辐射决定的
q_t = q_h-q_c_l-q_r_l; % 吸热管壁面总热流
% for ii = 1:floor(H*N_p/N_fp/deltaZ)
%     z = ii*deltaZ - H*floor(ii*deltaZ/H);
%     q_h(ii) = q_h(ii).*sin(z/H*pi); 
% end % 一种假设的分布：沿着吸热器周向均匀，轴向sin

% for ii = 1:floor(H*N_p/N_fp/deltaZ)
%     q_h(ii) = q_h(ii)* cos(floor(ii*deltaZ/H)/floor(N_p/N_fp)*pi/2);
% end % 一种假设的分布：沿着吸热器轴向均匀，周向cos

T_0 = ((epsilon_sky*T_sky^4 + epsilon_gr*T_amb^4)/(epsilon_gr+epsilon_sky))^0.25; % 环境温度

for z_num = 1:floor(H*N_p/N_fp/deltaZ)
   
    % z_num = 1
    for nn = 1:nnMax
        T_wall_old = T_wall(z_num,:); % 保存旧的壁面温度
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
        end % 最后一行，代表m=0

        Radb = zeros(Ns+1,1);
        for m=1:Ns
            Radb(m) = (kDelta(m,Ns+1)-F_(m,Ns+1))*T_Ns1(z_num)^4 + (kDelta(m,0)-F_(m,Ns+2))*T_0^4 ;
            for jj=1:Ns
                Radb(m) = (kDelta(m,jj)-F_(m,jj))*T_wall(z_num,jj)^4 + Radb(m);
            end
            Radb(m) = Radb(m) - F_(m,Ns+2).*q_h(z_num)/sig*alpha_;
        end

        Radb(Ns+1) = (kDelta(0,Ns+1)-F_(Ns+1,Ns+1))*T_Ns1(z_num)^4 + (kDelta(0,0)-F_(Ns+1,Ns+2))*T_0^4;
        for jj=1:Ns
            Radb(Ns+1) = (kDelta(0,jj)-F_(Ns+1,jj))*T_wall(z_num,jj)^4 + Radb(Ns+1);
        end
        Radb(Ns+1) = Radb(Ns+1) - F_(Ns+1,Ns+2)*q_h(z_num)/sig*alpha_;
        Radb = Radb*sig;
        q = RadA\Radb; % 求解辐射热流
        q = q/sig;
        q_0(z_num) = q(1); % 对应q_0
        q_j(z_num,:) = q(2:Ns+1)'; % 对应q_j


        q_c_l(z_num,:) = h.*(T_wall(z_num,:)-T_amb);
        % q_r_l(z_num,:) = -q_0(z_num) *2*B./(d_out.*deltaPhi.*(1:Ns)); % 吸热管壁面辐射热流 
        % q_r_l(z_num,:) = -q_0(z_num) *2*B./(d_out.*deltaPhi);
        q_r_l(z_num,:) = -q_0(z_num);
        q_t(z_num,:) = q_h(z_num) - q_c_l(z_num,:) - q_r_l(z_num,:);

        % if z_num == 1
        %     T_salt(z_num) = T_salt_0 + q_t(z_num,:)*(d_out/2*2*pi)*(deltaZ/2)/(m_*Cp_salt(T_salt_0)); % 计算盐温度
        % else
        %     T_salt(z_num) = T_salt(z_num-1) + q_t(z_num,:)*(d_out/2*2*pi)*(deltaZ)/(m_*Cp_salt(T_salt(z_num-1))); % 计算盐温度
        % end % 得到了熔融盐温度

        % 反算吸热管温度
        T_wall_mean(z_num) = mean(T_wall(z_num,:));
        if z_num == 1
            U = ((d_out*log(d_out/d_int))/(2*k_t) + R_foul*d_out/d_int + d_out/d_int/h_salt(T_salt_0,m_,d_int))^(-1);
            NTU = U*pi*d_out*(z_num-0.5)*(deltaZ*(z_num-0.5))/(m_*Cp_salt(T_salt_0));
        else
            U = ((d_out*log(d_out/d_int))/(2*k_t) + R_foul*d_out/d_int + d_out/d_int/h_salt(T_salt(z_num-1),m_,d_int))^(-1);
            NTU = U*pi*d_out*(z_num-0.5)*(deltaZ*(z_num-0.5))/(m_*Cp_salt(T_salt(z_num-1))); 
        end
        
        T_salt(z_num) = T_wall_mean(1) - (T_wall_mean(1)-T_salt_0)*exp(-NTU*(z_num-0.5)*deltaZ/(H*N_p/N_fp)); % 计算盐温度

        T_wall(z_num,:) = q_t(z_num,:)./U + T_salt(z_num); % 计算吸热管壁面温度

        % 计算保温层温度
        T_Ns1(z_num) = (kDelta(Ns+1,0)/eps0-(1/eps0 - 1)*F_(Ns+1,Ns+1))/sig*q_0(z_num); 
        for jj=1:Ns
            T_Ns1(z_num) = T_Ns1(z_num) + (kDelta(Ns+1,jj)/epsilon_t-(1/epsilon_t - 1)*F_(m,Ns-jj+1))/sig * q_j(z_num,jj);
        end % 对应q_j
        T_Ns1(z_num) = T_Ns1(z_num) - (kDelta(Ns+1,0) - F_(Ns+1,Ns+1))*T_0^4; % 对应q_0
        for jj=1:Ns
            T_Ns1(z_num) = T_Ns1(z_num) - (kDelta(Ns+1,jj)-F_(Ns+1,Ns-jj+1))*T_wall(z_num,jj)^4;
        end
        T_Ns1(z_num) = T_Ns1(z_num) + F_(Ns+1,Ns+1)*q_h(z_num)/sig*alpha_; % 对应q_h
        T_Ns1(z_num) = T_Ns1(z_num)^(0.25); % 保温层温度

        if max(abs(T_wall(z_num,:) - T_wall_old))< TOL
            fprintf('迭代完成，迭代次数为 %d, z_num = %d\n', nn, z_num);
            break; 
        elseif nn == nnMax
            fprintf('迭代未收敛，迭代次数为 %d, z_num = %d\n', nn, z_num);
            break;
        end
    end
end

save cauResult.mat
    






















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

function val = Cp_salt(T)
    val = 1443 + 0.172*(T-273.15);
end

function val = k_salt(T)
    val = 0.443 + 1.9e-4*(T-273.15); % W/m/K
end

function val = h_salt(T,m_,d_int)
    t = T - 273.15; % 摄氏度
    nu = 22.714 - 0.12*t + 2.281e-4*t^2 - 1.474e-7*t^3; % 动力粘度
    nu = nu*1e-3;
    Re = 4*m_/(pi*nu*d_int); % 雷诺数
    Pr = Cp_salt(T)*nu/k_salt(T);
    Nu = 0.023*Re^0.8*Pr^0.4;
    val = Nu*k_salt(T)/d_int; % 强制对流换热系数
end