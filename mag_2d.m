%2D magnetic modelling
%Mohammad Rheza Zamani
%Reference : Lowrie, W. 2007. Fundamental of Geophysics. Cambridge University Press, p. 330-332.
clear all;
clc;
%Input Paramater
mu0 = 1.26*10^-6; %H/m (Magnetic Permeability)
conts = (mu0)./(2*pi);
lengthx = 1000;
lengthz = 1000;
%Block dimenssion
dx = 20;
dz = 20;
%Block middle point 
dmx = dx/2;
dmz = dz/2;
%Number of vertical and horizontal block
nx = lengthx/dx;
nz = lengthz/dz;
%Total number of block
nb = nx*nz;
%Contrast suscepbility
V = zeros(nz,nx);
%V(24:25,1:12) = 0.9*10^-3;
%V(26:27,12:22) = 1*10^-3;
%V(28:29,22:32) = 0.7*10^-3;
%V(30:31,32:42) = 0.6*10^-3;
%V(32:33,42:50) = 0.5*10^-3;
V(1:12,24:25) = 1*10^-3;
V(12:22,26:27) = 0.8*10^-3;
%V(1:12,1:12) = 0.5*10^-3;

%Make block model
for i = 1 : nx
    x(i) = dx*i - dmx;
end
xx = repmat(x,nz,1);
for j = 1 : nz
    z(j) = dz*j;
end
z1 = [0 ;z'];
zz = repmat(z1,1,nx+1);
nb1 = length(zz)*nx;

%Kernell
for i = 1:nx
    for j = 1 : nb
        for k = 1 : nb1
            alpha1 = atan(((x(i)-xx(j))+dmx)./zz(k));
            alpha2 = atan(((x(i)-xx(j))-dmx)./zz(k));
            alpha3 = atan(((x(i)-xx(j))+dmx)./zz(k+1));
            alpha4= atan(((x(i)-xx(j))-dmx)./zz(k+1));
            Kernell(i,j) = conts.*((alpha1-alpha2)-(alpha3-alpha4));
        end
    end
end
%Calculated magnetic response
V_rs = reshape(V,nb,1);
dBdz = Kernell*V_rs;

%Plot kernel matrix
figure(1)
imagesc(Kernell)
set(gcf, 'Position', get(0, 'Screensize'));
ylabel('Observation Points','FontWeight','bold','FontSize',10)
xlabel('The Blocks','FontWeight','bold','FontSize',10)
cb = colorbar;
cb.Label.String = 'Value';
cb.Location = 'southoutside';
colormap(jet)

figure(2)
subplot(2,1,1)
plot(x,dBdz,'*-b')
xlabel('Distance (m)','FontWeight','bold','FontSize',10)
ylabel('\DeltaB_{z} (nT)','FontWeight','bold','FontSize',10)
title('Geomagnetic Response','FontWeight','bold','FontSize',10)
grid on
subplot(2,1,2)
s = pcolor(x,z,V);
s.FaceColor = 'interp';
xlabel('Distance (m)','FontWeight','bold','FontSize',10)
ylabel('Depth(m)','FontWeight','bold','FontSize',10)
title('Subsurface Model','FontWeight','bold','FontSize',10)
cb = colorbar;
cb.Label.String = 'DeltaM_{z} ';
cb.Location = 'southoutside';
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'Ydir','reverse')
colormap(jet)

data = [x' dBdz];
saveas(figure(1),'Kernell magnetic.png')
saveas(figure(2),'Model Magnetic.png')
writematrix(data,'Data forward magnetic.dat','Delimiter','tab')
