% Example where the Rosenbrock function is approximated
close all

% Rosenbrock function
rosenbrock = @(x,y) (1-x).^2 + 100*(y-x.^2).^2;

% Coarse grid for with sample points
N = 5;
x = linspace(-2,2,N);
y = linspace(-1,3,N);
[X,Y] = meshgrid(x,y);
Z = rosenbrock(X,Y);

% Fine grid for plotting and evaluation of errors
Nd = 20*N;
xd = linspace(-2,2,Nd);
yd = linspace(-1,3,Nd);
[Xd,Yd] = meshgrid(xd,yd);
Zd = rosenbrock(Xd,Yd);

% Plot function and sample points
figure
surf(Xd,Yd,Zd, 'EdgeColor','none','LineStyle','none')
hold on
stem3(X,Y,Z, 'color', 'black', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'MarkerFaceColor','black')
hold off
zlim([0,3000]);
view(210, 30);

% Sample function
d = DataTable;
for xi = x
   for yi = y
      d.add_sample([xi yi], rosenbrock(xi,yi));
   end
end

% Build approximations
approximator1 = BSpline(d, BSplineType.Linear);
approximator3 = BSpline(d, BSplineType.Cubic);

% Evaluate approximations and compute errors
approx1 = zeros(Nd,Nd);
error1 = approx1;

approx3 = zeros(Nd,Nd);
error3 = approx3;

i=1;
for xi = xd
    j = 1;
    for yi = yd
        exact = rosenbrock(xi,yi);
        
        approx1(i,j) = approximator1.eval([xi yi]);
        error1(i,j) = (approx1(i,j) - exact);
        approx3(i,j) = approximator3.eval([xi yi]);
        error3(i,j) = (approx3(i,j) - exact);
        
        j = j+1;
    end
   i = i+1;
end

% Plot approximations
figure
surf(xd, yd, approx1', 'EdgeColor','none','LineStyle','none')
zlim([0,3000]);
view(210, 30);

figure
surf(xd, yd, approx3', 'EdgeColor','none','LineStyle','none')
zlim([0,3000]);
view(210, 30);

% Plot errors
figure
surf(xd, yd, abs(error1)', 'EdgeColor','none','LineStyle','none')
%zlim([0,3000]);
view(210, 30);

figure
surf(xd, yd, abs(error3)', 'EdgeColor','none','LineStyle','none')
%zlim([0,3000]);
view(210, 30);

% Compute absolute errors
disp('Max error with linear spline:');
abserror1 = max(max(abs(error1)));
abserror1

disp('Max error with cubic spline:');
abserror3 = max(max(abs(error3)));
abserror3

% Compute relative errors
rangef = abs(max(max(Zd)) - min(min(Zd)));

disp('Max relative error with linear spline:');
relerror1 = abserror1/rangef;
relerror1

disp('Max relative error with cubic spline:');
relerror3 = max(max(abs(error3)))/rangef;
relerror3
