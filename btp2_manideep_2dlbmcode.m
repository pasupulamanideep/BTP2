%Compressible Lattice Boltzmann Method code for 2D Shock Problem
%BY PASUPULA MANIDEEP - 18AE30011
%For the coursework of PROJECT -II - AE47006 (BTP-II)

clc;
clear;
gamma = 1.4;
xlen = 104;
ylen = xlen;
x = 3:xlen-2;
y = 3:ylen-2;
D = 2;
b = 4;
j = 8;
k = 4;
rho = zeros(xlen,ylen);
p = zeros(xlen,ylen);
e = zeros(xlen,ylen);
u = zeros(xlen,ylen);
v = zeros(xlen,ylen);
phi = zeros(xlen,ylen);
c1_prime = zeros(xlen,ylen);
c2_prime = zeros(xlen,ylen);
d1 = zeros(xlen,ylen);
d2 = zeros(xlen,ylen);
c_bar_prime = zeros(xlen,ylen);
cj_prime_x = zeros(xlen,ylen,j);
cj_prime_y = zeros(xlen,ylen,j);
uk_prime = zeros(xlen,ylen,k);
vk_prime = zeros(xlen,ylen,k);
xcg = zeros(xlen);
ycg = zeros(ylen);
chi_jk = zeros(xlen,ylen,j,k);
m_jk = zeros(xlen,ylen,j,k);
xi_jk_x = zeros(xlen,ylen,j,k);
xi_jk_y = zeros(xlen,ylen,j,k);
zeta_jk = zeros(xlen,ylen,j,k);
d_jk = zeros(xlen,ylen,j,k);
rhoF = zeros(xlen,ylen,j,k);
rho_u_F = zeros(xlen,ylen,j,k);
rho_v_F = zeros(xlen,ylen,j,k);
rho_E_F = zeros(xlen,ylen,j,k);
rhoS = zeros(xlen,ylen,k);
rhok = zeros(xlen,ylen,k);
alphak = zeros(xlen,ylen,k);

%Initialising
[rho,p,u,v,e,phi,c1_prime,c2_prime,d1,d2,c_bar_prime] = initialising(xlen,ylen,gamma,D,b,rho,p,u,v,e,phi,c1_prime,c2_prime,d1,d2,c_bar_prime);

[xcg,ycg] = cgloc_calc(xcg,ycg,u,v,xlen,ylen);
[uk_prime,vk_prime] = uk_vk_prime_calc(uk_prime,vk_prime,xcg,ycg,xlen,ylen,u,v,k);
[cj_prime_x,cj_prime_y] = cj_prime_calc(cj_prime_x,cj_prime_y,c1_prime,c2_prime,xlen,ylen);
[rhok] = rhok_calc(rhok,rho,xlen,ylen,uk_prime,vk_prime);
[alphak] = alphak_calc(alphak,rhok,rho,xlen,ylen);
[chi_jk] = chijk_calc(chi_jk,D,cj_prime_x ,cj_prime_y,uk_prime,vk_prime,xlen,ylen,j,k);
[m_jk] = mjk_calc(chi_jk);
[xi_jk_x,xi_jk_y] = xijk_calc(xi_jk_x,xi_jk_y,chi_jk,cj_prime_x,cj_prime_y,u,v,xlen,ylen,j,k);
[zeta_jk] = zetajk_calc(zeta_jk, u,v,cj_prime_x,cj_prime_y,c_bar_prime,phi,chi_jk,j,k,xlen,ylen);
[d_jk] = d_jk_calc(d_jk,d1,d2,alphak,k,xlen,ylen);
[rhoF, rho_u_F, rho_v_F, rho_E_F] = F_values(rhoF, rho_u_F, rho_v_F, rho_E_F,m_jk,chi_jk,xi_jk_x,xi_jk_y,d_jk,xlen,ylen,j,k);


%streaming
[rhoS] = rhoS_calc(xlen,ylen,k,rhoF);
[rho_u_S,rho_v_S] = rho_u_v_S_calc(xlen,ylen,j,k,rho_u_F,rho_v_F);
[rho_E_S] = rho_E_S_calc(xlen,ylen,k,rho_E_F);

[rho_new,rho_u_new,rho_v_new,rho_E_new,u_new,v_new,E_new,e_new,phi_new,c1_prime_new,c2_prime_new,d1_new,d2_new,c_bar_prime_new,p_new,xcg_new,ycg_new] = new_value_calc(rhoS,rho_u_S,rho_v_S,rho_E_S,k,xlen,ylen,D,b,gamma);

xcg = xcg_new;
ycg = ycg_new;
rho = rho_new;
p=p_new;
u = u_new;
v = v_new;
phi = phi_new;
c1_prime = c1_prime_new;
c2_prime = c2_prime_new;
d1 = d1_new;
d2 = d2_new;
c_bar_prime = c_bar_prime_new;


figure
hh = contourf(1:xlen,1:ylen,rho_new,10);
hold on
xlabel('x')
ylabel('y')
title('Plot of Density contours')

%keeping number of iterations as 100

for ii=1:100
[uk_prime,vk_prime] = uk_vk_prime_calc(uk_prime,vk_prime,xcg,ycg,xlen,ylen,u,v,k);
[cj_prime_x,cj_prime_y] = cj_prime_calc(cj_prime_x,cj_prime_y,c1_prime,c2_prime,xlen,ylen);
[rhok] = rhok_calc(rhok,rho,xlen,ylen,uk_prime,vk_prime);
[alphak] = alphak_calc(alphak,rhok,rho,xlen,ylen);
[chi_jk] = chijk_calc(chi_jk,D,cj_prime_x ,cj_prime_y,uk_prime,vk_prime,xlen,ylen,j,k);
[m_jk] = mjk_calc(chi_jk);
[xi_jk_x,xi_jk_y] = xijk_calc(xi_jk_x,xi_jk_y,chi_jk,cj_prime_x,cj_prime_y,u,v,xlen,ylen,j,k);
[zeta_jk] = zetajk_calc(zeta_jk, u,v,cj_prime_x,cj_prime_y,c_bar_prime,phi,chi_jk,j,k,xlen,ylen);
[d_jk] = d_jk_calc(d_jk,d1,d2,alphak,k,xlen,ylen);
[rhoF, rho_u_F, rho_v_F, rho_E_F] = F_values(rhoF, rho_u_F, rho_v_F, rho_E_F,m_jk,chi_jk,xi_jk_x,xi_jk_y,d_jk,xlen,ylen,j,k);
%streaming
[rhoS] = rhoS_calc(xlen,ylen,k,rhoF);
[rho_u_S,rho_v_S] = rho_u_v_S_calc(xlen,ylen,j,k,rho_u_F,rho_v_F);
[rho_E_S] = rho_E_S_calc(xlen,ylen,k,rho_E_F);
[rho_new,rho_u_new,rho_v_new,rho_E_new,u_new,v_new,E_new,e_new,phi_new,c1_prime_new,c2_prime_new,d1_new,d2_new,c_bar_prime_new,p_new] = new_value_calc(rhoS,rho_u_S,rho_v_S,rho_E_S,k,xlen,ylen,D,b,gamma);

xcg = xcg_new;
ycg = ycg_new;
rho = rho_new;
p = p_new;
u = u_new;
v = v_new;
phi = phi_new;
c1_prime = c1_prime_new;
c2_prime = c2_prime_new;
d1 = d1_new;
d2 = d2_new;
c_bar_prime = c_bar_prime_new;

hh = contourf(1:xlen,1:ylen,rho_new,10);
pause(.1)
end
hold off


%%%%%%%%%%Checking the properties that should be satisfied by equlibrium
%%%%%%%%%%distribution - equation 26 & 27 in sun-2003 paper
uk = zeros(xlen,ylen);
vk = zeros(xlen,ylen);
for tt = 1:4
uk(x,y,tt) = u(x,y) + uk_prime(x,y,tt);
vk(x,y,tt) = v(x,y) + vk_prime(x,y,tt);  
end

%condition1 check
rhouk = zeros(xlen,ylen);
rhouk_diff = zeros(xlen,ylen);
for tt = 1:4
    rhouk(x,y) = rhouk(x,y) + rhok(x,y,tt).*uk(x,y,tt);
end
rhouk_diff(x,y) = rhouk(x,y) - rho(x,y).*u(x,y);

%condition2 check
rhovk = zeros(xlen,ylen);
rhovk_diff = zeros(xlen,ylen);
for tt = 1:4
    rhovk(x,y) = rhovk(x,y) + rhok(x,y,tt).*vk(x,y,tt);
end
rhovk_diff(x,y) = rhovk(x,y) - rho(x,y).*v(x,y);

%condition3 check
rhok_sum = zeros(xlen,ylen);
rhok_sum_diff = zeros(xlen,ylen);
for tt = 1:4
    rhok_sum(x,y) = rhok_sum(x,y) + rhok(x,y,tt);
end
rhok_sum_diff(x,y) = rhok_sum(x,y) - rho(x,y);

%condition4 check
rhovk_prime = zeros(xlen,ylen);
for tt = 1:4
    rhovk_prime(x,y) = rhovk_prime(x,y) + rhok(x,y,tt).*vk_prime(x,y,tt);
end

%condition5 check
rhouk_prime = zeros(xlen,ylen);
for tt = 1:4
    rhouk_prime(x,y) = rhouk_prime(x,y) + rhok(x,y,tt).*uk_prime(x,y,tt);
end

%condition6 check
alphak_sum = zeros(xlen,ylen);
for tt = 1:4
    alphak_sum(x,y) = alphak_sum(x,y) + alphak(x,y,tt);
end




%%%%Functions

function [rho,p,u,v,e,phi,c1_prime,c2_prime,d1,d2,c_bar_prime] = initialising(xlen,ylen,gamma,D,b,rho,p,u,v,e,phi,c1_prime,c2_prime,d1,d2,c_bar_prime)
xlen_f = xlen - 4;
ylen_f = ylen - 4;
for m = 3:xlen-2
    for n = 3:ylen-2
        if (m > (0.25*xlen_f +2) && m < (0.75*xlen_f +2)) && (n > (0.25*ylen_f +2) && n < (0.75*ylen_f +2)) %inner box
            rho(m,n) = 0.125;
            p(m,n) = 0.025;
            u(m,n) = 0;
            v(m,n) = 0;
        else %outer box
            rho(m,n) = 1;
            p(m,n) = 0.25;
            u(m,n) = 0;
            v(m,n) = 0;            
        end
    end
end

m = 3:xlen-2;
n = 3:ylen-2;

e(m,n) = p(m,n)./((gamma-1).*rho(m,n));
phi(m,n) = (1-(gamma-1)*D/2).*e(m,n);
c1_prime(m,n) = floor(sqrt(((gamma-1)*D).*e(m,n)));
c2_prime(m,n) = c1_prime(m,n)+1;
d1(m,n) = rho(m,n).*((c2_prime(m,n).^2 - ((gamma-1)*D).*e(m,n))./(b.*(c2_prime(m,n).^2 - c1_prime(m,n).^2)));
d2(m,n) = rho(m,n).*((c1_prime(m,n).^2 - ((gamma-1)*D).*e(m,n))./(b.*(c1_prime(m,n).^2 - c2_prime(m,n).^2)));
c_bar_prime(m,n) = sqrt(b.*(d1(m,n).* c1_prime(m,n).^2 + d2(m,n).* c2_prime(m,n).^2)./rho(m,n));
end



function [xcg,ycg] = cgloc_calc(xcg,ycg,u,v,xlen,ylen)
for x = 3:xlen-2
for y = 3:ylen-2

xcg(x,y) = x + u(x,y);
ycg(x,y) = y + v(x,y);

end
end
end



function [uk_prime,vk_prime] = uk_vk_prime_calc(uk_prime,vk_prime,xcg,ycg,xlen,ylen,~,~,~)
for x=3:xlen-2
    for y = 3:ylen-2
        uk_prime(x,y,1) = -(xcg(x,y) - x);
        vk_prime(x,y,1) = -(ycg(x,y) - y);
        uk_prime(x,y,4) = uk_prime(x,y,1);
        uk_prime(x,y,3) = 1 - abs(uk_prime(x,y,1));
        uk_prime(x,y,2) = uk_prime(x,y,3);
        vk_prime(x,y,2) = vk_prime(x,y,1);
        vk_prime(x,y,3) = 1 - abs(vk_prime(x,y,2));
        vk_prime(x,y,4) = vk_prime(x,y,3);
    end
end
end



function [cj_prime_u,cj_prime_v] = cj_prime_calc(cj_prime_u,cj_prime_v,c1_prime,c2_prime,xlen,ylen)
m = 3:xlen-2;
n = 3:ylen-2;
dir_u = [1 0 1 0 1 0 1 0]';
for j=1:8
    if j<=4
        cj_prime_u(m,n,j) = c1_prime(m,n).*dir_u(j);
    else 
        cj_prime_u(m,n,j) = c2_prime(m,n).*dir_u(j);
    end
end
dir_v = [0 1 0 1 0 1 0 1]';
for j=1:8
    if j<=4
        cj_prime_v(m,n,j) = c1_prime(m,n).*dir_v(j);
    else 
        cj_prime_v(m,n,j) = c2_prime(m,n).*dir_v(j);
    end
end
end



function [rhok] = rhok_calc(rhok,rho,xlen,ylen,uk_prime,vk_prime)
m = 3:xlen-2;
n = 3:ylen-2;
    rhok(m,n,1)=rho(m,n).*abs(uk_prime(m,n,3)).*abs(vk_prime(m,n,3));
    rhok(m,n,2)=rho(m,n).*abs(uk_prime(m,n,4)).*abs(vk_prime(m,n,4));
    rhok(m,n,3)=rho(m,n).*abs(uk_prime(m,n,1)).*abs(vk_prime(m,n,1));
    rhok(m,n,4)=rho(m,n).*abs(uk_prime(m,n,2)).*abs(vk_prime(m,n,2));
end



function [alphak] = alphak_calc(alphak,rhok,rho,xlen,ylen)
m = 3:xlen-2;
n = 3:ylen-2;
for k=1:4
    alphak(m,n,k) = rhok(m,n,k)./rho(m,n);
end
end



function [chijk] = chijk_calc(chijk,~,cj_prime_u,cj_prime_v,uk_prime,vk_prime,xlen,ylen,j,k)
m = 3:xlen-2;
n = 3:ylen-2;


for jj = 1:j
    for kk = 1:k
          chijk(m,n,jj,kk) = (cj_prime_u(m,n,jj).*uk_prime(m,n,kk) + cj_prime_v(m,n,jj).*vk_prime(m,n,kk))./(cj_prime_u(m,n,jj).^2 + cj_prime_v(m,n,jj).^2 + 10^-6);
    end
end
end



function [m_jk] = mjk_calc(chi_jk)
m_jk = 1 - chi_jk;
end



function [xi_jk_u,xi_jk_v] = xijk_calc(xi_jk_u,xi_jk_v,chi_jk,cj_prime_u,cj_prime_v,u,v,xlen,ylen,j,k)
m = 3:xlen-2;
n = 3:ylen-2;

for kk=1:k
    for jj = 1:j
        xi_jk_u(m,n,jj,kk) = u(m,n) + cj_prime_u(m,n,jj) - chi_jk(m,n,jj,kk).*u(m,n);
        xi_jk_v(m,n,jj,kk) = v(m,n) + cj_prime_v(m,n,jj) - chi_jk(m,n,jj,kk).*v(m,n);
    end
end
end



function [d_jk] = d_jk_calc(d_jk,d1,d2,alphak,k,xlen,ylen)
for x=3:xlen-2
    for y=3:ylen-2
for m=1:4
    for n=1:k
        d_jk(x,y,m,n) = alphak(x,y,n)*d1(x,y);
    end
end
for m=5:8
    for n=1:k
        d_jk(x,y,m,n) = alphak(x,y,n)*d2(x,y);
    end
end
    end
end
end   



function [zeta_jk] = zetajk_calc(zeta_jk, u,v,cj_prime_u,cj_prime_v,c_bar_prime,phi,~,j,k,xlen,ylen)
m = 1:xlen;
n = 1:ylen;

for kk = 1:k
    for jj = 1:j
% zeta_jk(m,n,jj,kk) = (1/2).*((u(m,n).^2 + v(m,n).^2) + 2.*(cj_prime_u(m,n,jj).*u(m,n) + cj_prime_v(m,n,jj).*v(m,n)) + c_bar_prime(m,n).^2)- chi_jk(m,n,jj,kk).*((1/2).*((u(m,n).^2 + v(m,n).^2) + c_bar_prime(m,n).^2)) + phi.*(1-chi_jk(m,n,jj,kk));           
    zeta_jk(m,n,jj,kk) = (1/2).*((u(m,n).^2 + v(m,n).^2) + 2.*(cj_prime_u(m,n,jj).*u(m,n) + cj_prime_v(m,n,jj).*v(m,n)) + c_bar_prime(m,n).^2) + phi; 
    end
end
end



function [rhoF, rho_u_F, rho_v_F, rho_E_F] = F_values(rhoF, rho_u_F, rho_v_F, rho_E_F,m_jk,chi_jk,xi_jk_x,xi_jk_y,d_jk,xlen,ylen,j,k)
for x=3:xlen-2
    for y=3:ylen-2
        for m=1:j
            for n=1:k
                rhoF(x,y,m,n) = d_jk(x,y,m,n)*m_jk(x,y,m,n);
                rho_u_F(x,y,m,n) = d_jk(x,y,m,n)*xi_jk_x(x,y,m,n);
                rho_v_F(x,y,m,n) = d_jk(x,y,m,n)*xi_jk_y(x,y,m,n);
                rho_E_F(x,y,m,n) = d_jk(x,y,m,n)*chi_jk(x,y,m,n);
            end
        end
    end
end
end


function [rhoS] = rhoS_calc(xlen,ylen,k,rhoF)
rhoS = zeros(xlen,ylen,k);
for x=3:xlen-2
    for y=3:ylen-2
        for n=1:k
        rhoS(x,y,n) = rhoF(x+1,y,3,n)+rhoF(x+2,y,7,n)+rhoF(x-1,y,1,n)+rhoF(x-2,y,5,n)+rhoF(x,y+1,2,n)+rhoF(x,y+2,6,n)+rhoF(x,y-1,4,n)+rhoF(x,y-2,8,n);
        end
       
    end
end
end


function [rho_u_S,rho_v_S] = rho_u_v_S_calc(xlen,ylen,~,k,rho_u_F,rho_v_F)
rho_u_S = zeros(xlen,ylen,k);
rho_v_S = zeros(xlen,ylen,k);
for x=3:xlen-2
    for y=3:ylen-2
        for n=1:k
        rho_u_S(x,y,n) = (rho_u_F(x+1,y,3,n))+rho_u_F(x+2,y,7,n)+rho_u_F(x-1,y,1,n)+rho_u_F(x-2,y,5,n);
        rho_v_S(x,y,n) = (rho_v_F(x,y+1,2,n))+rho_v_F(x,y+2,6,n)+rho_v_F(x,y-1,4,n)+rho_v_F(x,y-2,8,n);
        end
        
    end
end
end



function [rho_E_S] = rho_E_S_calc(xlen,ylen,k,rho_E_F)
rho_E_S = zeros(xlen,ylen,k);
for x=3:xlen-2
    for y=3:ylen-2
        for n=1:k
            rho_E_S(x,y,n) = rho_E_F(x+1,y,3,n)+rho_E_F(x+2,y,7,n)+rho_E_F(x-1,y,1,n)+rho_E_F(x-2,y,5,n)+rho_E_F(x,y+1,2,n)+rho_E_F(x,y+2,6,n)+rho_E_F(x,y-1,4,n)+rho_E_F(x,y-2,8,n);
        end
    end
end
end



function [rho_new,rho_u_new,rho_v_new,rho_E_new,u_new,v_new,E_new,e_new,phi_new,c1_prime_new,c2_prime_new,d1_new,d2_new,c_bar_prime_new,p_new,xcg_new,ycg_new] = new_value_calc(rhoS,rho_u_S,rho_v_S,rho_E_S,k,xlen,ylen,D,b,gamma)
rho_new = zeros(xlen,ylen);
rho_u_new = zeros(xlen,ylen);
rho_v_new = zeros(xlen,ylen);
rho_E_new = zeros(xlen,ylen);
c1_prime_new = zeros(xlen,ylen);
c2_prime_new = zeros(xlen,ylen);
d1_new = zeros(xlen,ylen);
d2_new = zeros(xlen,ylen);
c_bar_prime_new = zeros(xlen,ylen);
phi_new = zeros(xlen,ylen);
p_new = zeros(xlen,ylen);
e_new = zeros(xlen,ylen);
xcg_new = zeros(xlen,ylen);
ycg_new = zeros(xlen,ylen);
for x=3:xlen-2
    for y=3:ylen-2
        for nn=1:k
            rho_new(x,y) = rho_new(x,y)+rhoS(x,y,nn);
            rho_u_new(x,y) = rho_u_new(x,y)+ rho_u_S(x,y,k);
            rho_v_new(x,y) = rho_v_new(x,y)+ rho_v_S(x,y,k);
            rho_E_new(x,y) = rho_E_new(x,y)+ rho_E_S(x,y,k);
        end
    end
end
        u_new = rho_u_new./rho_new;
        v_new = rho_v_new./rho_new;
        E_new = rho_E_new./rho_new;

for x = 3:xlen-2
for y = 3:ylen-2

xcg_new(x,y) = x + u_new(x,y);
ycg_new(x,y) = y + v_new(x,y);

end
end
        
        
        
%other values
mm = 3:xlen-2;
nn = 3:ylen-2;

e_new(mm,nn) = E_new(mm,nn)-(0.5).*(u_new(mm,nn).^2 + v_new(mm,nn).^2);
p_new(mm,nn) = (gamma - 1).*rho_new(mm,nn).*e_new(mm,nn);
phi_new(mm,nn) = (1-(gamma-1)*D/2).*e_new(mm,nn);
c1_prime_new(mm,nn) = floor(sqrt(((gamma-1)*D).*e_new(mm,nn)));
c2_prime_new(mm,nn) = c1_prime_new(mm,nn)+1;
d1_new(mm,nn) = rho_new(mm,nn).*((c2_prime_new(mm,nn).^2 - ((gamma-1)*D).*e_new(mm,nn))./(b.*(c2_prime_new(mm,nn).^2 - c1_prime_new(mm,nn).^2)));
d2_new(mm,nn) = rho_new(mm,nn).*((((gamma-1)*D).*e_new(mm,nn)- c1_prime_new(mm,nn).^2)./(b.*(c2_prime_new(mm,nn).^2 - c1_prime_new(mm,nn).^2)));
c_bar_prime_new(mm,nn) = sqrt(b.*(d1_new(mm,nn).* c1_prime_new(mm,nn).^2 + d2_new(mm,nn).* c2_prime_new(mm,nn).^2)./rho_new(mm,nn));
end

