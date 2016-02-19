% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

% B-spline regression of noisy sine function
close all

% Noisy sine function
myfunc = @(x) sin(x) + normrnd(0, 0.5);

% Coarse grid for with sample points
N = 1000;
xstart = 0;
xend = 4*pi;
x = linspace(xstart,xend,N);
y = zeros(1,N);
for i = 1:N
    y(i) = myfunc(x(i));
end

% Fine grid for plotting and evaluation of errors
Nd = N;
xd = linspace(xstart,xend,Nd);
yd = zeros(Nd,1);
for i = 1:Nd
    yd(i) = myfunc(xd(i));
end

% Sample function
d = DataTable;
for i = 1:N
  d.add_sample(x(i), y(i));
end

% Build approximations
approx1 = BSpline(d, BSplineType.Linear);
approx2 = BSpline(d, BSplineType.Quadratic);
approx3 = BSpline(d, BSplineType.Cubic);

% Evaluate approximations
yad1 = zeros(Nd,1);
yad2 = zeros(Nd,1);
yad3 = zeros(Nd,1);
% error1 = approx1;

for i = 1:Nd
    yad1(i) = approx1.eval(xd(i));
    yad2(i) = approx2.eval(xd(i));
    yad3(i) = approx3.eval(xd(i));
end

% Plot sample points and approximations
figure
plot(x,y,'.k')
hold on
plot(xd,yad1,'--k')
plot(xd,yad2,'-.k')
plot(xd,yad3,'-k')
hold off
legend('Sine', 'Linear', 'Quadratic', 'Cubic');
