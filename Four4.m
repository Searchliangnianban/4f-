close all;
clear all;
clc;
L1=0.0005; %side length 
M=1000; %number of samples 
dx1=L1/M; %src sample interval 
x1=-L1/2:dx1:L1/2-dx1; %src coords 
y1=x1;

lambda=0.5*10^-6; %wavelength 
k=2*pi/lambda; %wavenumber 
w1=0.000022; %source half width (m) 
w2=0.000011;
z=0.5*10^-3;
% z0=0.5*10^-3;
P=1e-5*0.5;
%sita=60;
b = lambda*z/P;

[X1,Y1]=meshgrid(x1,y1); 
R = sqrt(X1.^2+Y1.^2);


%u1 = circ(R/w1);%物面

% u00=ifftshift(ifft2(fftshift(u1)));

%u1=rect(X1/(2*w1)).*rect(Y1/(2*w1));
%u1=1/2*(1+cos(2*pi*X1/(P))).*rect(X1/(w1)).*rect(Y1/(w1));
%u1=1/2*(1+cos(2*pi*X1/(P)+pi/2));
%u1 = sin(X1) + cos(Y1);
u1=rect((X1-b)/(2*w1)).*rect(Y1/(2*w2))+rect((X1+b)/(2*w2)).*rect(Y1/(2*w1));


lens1 = exp(-1i*k*(X1.^2+Y1.^2)/(2*z));
lens2 = exp(-1i*k*(X1.^2+Y1.^2)/(2*z));

I1=abs(u1.^2);

H=1/2*(1+cos(2*pi*X1/(P)));

u21 = propIR(u1,L1,lambda,z);
u22 = u21.*lens1;
u3 = propIR(u22,L1,lambda,z).*H;%频谱面
%[X2,Y2]=meshgrid(x2,x2);

I2 = abs(u3.^2);%频谱面强度图



u41 = propIR(u3,L1,lambda,z);
u42 = u41.*lens2;
u5 = propIR(u42,L1,lambda,z);

I3 = abs(u5.^2);%像平面
I4 = abs(u41.^2);

figure(1) 
imagesc(x1,y1,I1); 
%surf(x1,y1,abs(u1.^2)); 
axis square; axis xy; 
colorbar();
colormap("hot"); xlabel('x (m)'); ylabel('y (m)'); 
title('物平面');

figure(2) 
imagesc(x1,y1,I2);
%surf(x2,y2,abs(u2.^2)); 
axis square; axis xy; 
colorbar();
colormap("hot"); xlabel('x (m)'); ylabel('y (m)'); 
title('空间频率频谱面');

figure(3) 
imagesc(x1,y1,I3);
%surf(x3,y3,abs(u3.^2)); 
axis square; axis xy; 
colorbar();
colormap("hot"); xlabel('x (m)'); ylabel('y (m)'); 
title('像面');

figure(4) 
mesh(x1,y1,I2); 
axis square; axis xy; 
colorbar();
colormap("hot"); xlabel('x (m)'); ylabel('y (m)');zlabel('z'); 
title('3D频谱面');

% figure(5) %irradiance profile 
% plot(x1,I2(M/2+1,:)); 
% xlabel('x (m)'); ylabel('Irradiance'); 
% title(['z= ',num2str(z),' m']);


% figure(6) 
% imagesc(x1,y1,I4);
% %surf(x3,y3,abs(u3.^2)); 
% axis square; axis xy; 
% colorbar();
% colormap("hot"); xlabel('x (m)'); ylabel('y (m)'); 
% title('U21面');