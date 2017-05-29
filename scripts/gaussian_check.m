% random gaussian check
clc;clear all;close all

x0=0;
y0=0;

rad = 1;
A = 1;
% x = -5:0.01:5;
% y = x;
[x,y] = meshgrid(-5:0.01:5,-5:0.01:5);
gauss = A *exp( - (x.^2 + y.^2)/rad);
gauss(ceil(end/2)+1:end,:) = 0;
%%
% figure(); 
% surf(x,y,gauss)

%% 
% P = [0.1 0.3 0.4 0.2]
% X = [1 2 3 4]
% 
% C = cumsum(P);
% f = X(1+sum(C(end)*rand>C))
%%
C = cumsum(gauss(:))/sum(gauss(:));
X = x(:);
Y = y(:);
placeInd = 1+sum(C(end)*rand>C);
x_place = X(placeInd );
y_place = Y(placeInd );
% figure();plot(C)
% f = 

figure(); 
surf(x,y,gauss)
hold on
scatter3(x_place,y_place,1,'fill','r')
view([0,90])