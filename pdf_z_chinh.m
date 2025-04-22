 clc;clear;
% a=4;
% Fx = @(x) -2/a^3*x.^3+3/a^2*x.^2;
% fx = @(x) -6/a^3*x.^2+6/a^2*x;
R=10;
phi_1_2 = pi/3;
phi_0 = pi/4;
A_r = 5*10^(-5);
n=1.33;
FOV = pi/3;
T_s =1;
c_av = 0.0352;
H_1 = cos(phi_0)*R;
H_2 = R - H_1;
a = sin(phi_0)*R;



%Random gia tri z

% z =0:0.01:10;
 F_z_1 =  @(z) 3/2 * (1/(H_2*R^2))*((z.^3)/3*a^2/H_1^2);
 f_z_1 = @(z) 3/2 * (1/(H_2*R^2))*((z.^2)*a^2/H_1^2);
 F_z_2 = @(z) 3/2 * (1/(H_2*R^2))*(R^2*z-z.^3/3);
 f_z_2 = @(z) 3/2 * (1/(H_2*R^2))* (R^2 -z.^2);
% F_dist = makedist('PiecewiseLinear', 'z', z, 'F_z', F_z);
diem_1 = 3/2 * (1/(H_2*R^2))*((H_1^3)/3*a^2/H_1^2);
%d = 3/2 * (1/(H_2*R^2))*(R^2*H_1-H_1^3/3);

%N = 10000;
%x = zeros(1,N);
%u = zeros(1,N);
%y = linspace(0,H_1,30);
%y_1 = linspace(H_1,R,30);
%{
for i = 1:N
    u(i) = rand;
    if u(i) < diem_1
        x(i) = fzero(@(z) (F_z_1(z)-u(i)),H_1);
        if x(i)<0
            x(i )= NaN ;
         
        end
    end     
   if u(i)>diem_1
       h = u(i)-diem_1;
      x(i) = fzero(@(z)  (integral(f_z_2,H_1,z)-h),H_1);
        
       if x(i)<0
            x(i) = NaN;
         
        end
    end
    
    
    
    
   % x(i) = fzero(@(z) (F_z(z)-u(i)),1);
end
%}
    u = rand;
    if u < diem_1
        z_1 = fzero(@(z) (F_z_1(z)-u),H_1);
        if z_1<0
            z_1= NaN ;
         
        end
    end     
   if u > diem_1
       q = u-diem_1;
      z_1 = fzero(@(z)  (integral(f_z_2,H_1,z)-q),H_1);
        
       if z_1<0
            z_1 = NaN;
         
        end
    end
   



%plot(y,f_z_1(y))
%hold on
%{plot(y_1,f_z_2(y_1))
%hold on
%histogram(x,50,'Normalization','pdf')
%





%Random gia tri r
f_r = @(r) (3*r.*sqrt(R^2-r.^2)/(H_2*R^2))-(3*r.^2*H_1)/(H_2*R^2*a);
F_r = @(r) arrayfun(@(z) integral(f_r, 0, z), r);
u_1 = rand;
r_1 = fzero(@(r) (F_r(r)-u_1),2.5);

%Tinh goc phi
phi = atan(r_1/z_1);
%Cac gia tri khac
d = sqrt(z_1^2+r_1^2);
m = -log(2)/log(cos(phi_1_2));
A_eff = (A_r * cos(phi) * n^2)/(sin(FOV))^2;
%h = ((m+1)*(cos(phi))^m*A_eff*exp(-d*c_av))/(2*pi*d^2);
% Cac an lien quan den channel gain 
B = ((m+1)*A_r*n^2)/(2*pi*(sin(FOV))^2);
uu_1 = (4*m+6)/(m+1);
uu_2 = (3*m+5)/(m+1);
h_1 = (B*exp(-c_av*R))/R^2;
h_min = (B*exp(-c_av*R)*H_1^(m+1))/(R^(m+3));
e_1 = (3*(m+3)^uu_1)/(4*(m+1)*H_2*R^2*c_av^uu_2*B^(1/(m+1)));
uu_3 = (m+3)/(m+1);
e_2 = (R*c_av)/(m+3);
e_3 = (R^(m+3)/(B*exp(-R*c_av)))^(1/(m+1));
% SNR
y_t = 10^5;
y_1 = y_t *h_1^2;
y_min = y_t*h_min^2;
y_max = inf;

z_max = @(y) e_3 *sqrt(y./y_t).^(1/(m+1));
z_sao = @(y) (2/c_av)* lambertw(0, c_av/2*sqrt(B./sqrt(y./y_t)));
z_min = @(y)  2*H_1/(c_av*R)*lambertw(0,c_av/2*sqrt(B./sqrt(y./y_t)*(H_1/R)^(m+1)));
z_min_1 = @(y) z_min(sqrt(y)/sqrt(y_t));
denta_1 = @(y)   denta_tong(3,R*c_av*z_min(y)./(2*H_1),c_av*z_sao(y)./2,uu_3,m) + denta_tong(4,R*c_av*z_min(y)./(2*H_1),c_av*z_sao(y)./2,uu_3,m);
denta_2 = @(y)   denta_tong(3,R*c_av*z_min(y)./(2*H_1),e_2*(m+3)/2,uu_3,m) + denta_tong(4,R*c_av*z_min(y)./(2*H_1),e_2*(m+3)/2,uu_3,m);


f_y_1 = @(y) 3*(H_1^2*z_sao(y).^3 - R^2 * z_min(y).^3)./(4*H_1^2*H_2*R^2*(m+1)) .* y.^(-1) - (e_1*denta_1(y))./(2*y_t^(1/(2*m+2))*(m+1)).*y.^(-(2*m+1)/(2*m+2));

f_y_2 = @(y) -(3*z_min(y).^3./(4*H_1^2*H_2*(m+1))).*y.^-1 +(3*e_3/(4*H_2*y_t^(1/(2*m+2))*(m+1)) - e_1*denta_2(y)./(2*y_t^(1/(2*m+2))*(m+1))).*y.^(-(2*m+1)/(2*m+2));



F_y_1 = @(y) 1+ (H_1^2*z_sao(y).^3 - R^2*z_min(y).^3) ./ (2*H_1^2*H_2*R^2) - e_1*denta_1(y) ./ y_t^(1/(2+2*m)) .* y.^(1/(2+2*m));
F_y_2 = @(y) (-2*H_1^3 + 3*H_1^2*z_max(y) - z_min(y).^3)/(2*H_1^2*H_2) - e_1 * denta_2(y) ./ y_t^(1/(2+2*m)) .* y.^(1/(2+2*m));

f_h_1 = @(h) 3*(H_1^2*z_sao(h).^3 - R^2 * z_min(h).^3)./(2*H_1^2*H_2*R^2*(m+1)) .* h.^(-1) - (e_1*denta_1(h))./(m+1).*h.^(-(m)/(m+1));
f_h_2 = @(y) -(3*z_min(y).^3./(2*H_1^2*H_2*(m+1))).*y.^-1 +(3*e_3/(2*H_2*(m+1)) - e_1*denta_2(y)./(m+1)).*y.^(-(m)/(m+1));
F_h_1 = @(y) 1+ (H_1^2*z_sao(y).^3 - R^2*z_min(y).^3) ./ (2*H_1^2*H_2*R^2) - e_1*denta_1(y) .* y.^(1/(1+m));
F_h_2 = @(y) (-2*H_1^3 + 3*H_1^2*z_max(y) - z_min(y).^3)/(2*H_1^2*H_2) - e_1 * denta_2(y) .* y.^(1/(1+m));


%plot(y,f_y_1(y));

% plot(10*log10(y),f_y_2(y));
% hold on;
%}
% tich_phan = integral(f_y_2,y_min,y_1) + integral(f_y_1,y_1,y_max);
% 
% f_y_y_1 =@(y) f_y_1(y)./tich_phan;
% f_y_y_2 =@(y) f_y_2(y)./tich_phan;

%y = y_1:0.1:3200;
%plot(10*log10(y),f_y_y_1(y));
% hold on;


% Ve simulation


