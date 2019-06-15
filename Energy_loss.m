clear all; clc; close all;
 
% BASE DATA
filename = 'projectData3.xlsx';
file = xlsread(filename);
X = file(:, 1)';
Y = file(:, 2)';
p0 = plot(X, Y, '-o', 'MarkerSize', 5);
xlabel('Displacement(10^-3 in.)')
ylabel('Force(lb.)')
grid on
hold on
 
% LINEAR LEAST SQUARES
coefL = lspoly(X, Y, 1);
XX = 0:.01:max(X);
YY = zeros(1, length(XX));
for i = 1:length(XX)
    YY(i) = horner(coefL, XX(i));
end
p1 = plot(XX, YY);
 
% QUADRATIC LEAST SQUARES
coefQ = lspoly(X, Y, 2);
XXX = 0:.01:max(X);
YYY = zeros(1, length(XXX));
for i = 1:length(XXX)
    YYY(i) = horner(coefQ, XXX(i));
end
p2 = plot(XXX, YYY);
 
% CUBIC SPLINES
N = length(X);
n = N-1;
Unknowns = 4*n;
% Create the Y vector
B = zeros(4*N-4, 1);
% H0 is the equations of the zero derivatives (interpolation)
Num_row = 2*N - 2;
H0 = zeros(Num_row, Unknowns);
for idx = 1:(Num_row/2)
    col = (idx-1)*4 + 1;
    row = (idx-1)*2 + 1;
    for jdx = col:(col+3)
        H0(row, jdx)   = X(idx)^(3-jdx+col);
        H0(row+1, jdx) = X(idx+1)^(3-jdx+col);
    end
    B(row, 1)   = Y(idx);
    B(row+1, 1) = Y(idx+1);
end
% K1 is the equations of the first derivatives
K1 = zeros(N-2, Unknowns);
for idx = 2:(N-1)
    col = 1 + (idx-2)*4;
    K1(idx-1,col)   = 3*X(idx)^2;
    K1(idx-1,col+1) = 2*X(idx);
    K1(idx-1,col+2) = 1;
    K1(idx-1,col+3) = 0;
    K1(idx-1,col+4) = -3*X(idx)^2;
    K1(idx-1,col+5) = -2*X(idx);
    K1(idx-1,col+6) = -1;
    K1(idx-1,col+7) = 0;
end
% K2 is the equations of the second derivatives
K2 = 0*K1;
for idx = 2:(N-1)
    col = 1 + (idx-2)*4;
    K2(idx-1,col)   = 6*X(idx);
    K2(idx-1,col+1) = 2;
    K2(idx-1,col+2) = 0;
    K2(idx-1,col+3) = 0;
    K2(idx-1,col+4) = -6*X(idx);
    K2(idx-1,col+5) = -2;
    K2(idx-1,col+6) = 0;
    K2(idx-1,col+7) = 0;
end
% D2 is equations of the second derivative at the endpoints,
% which we set to zero.
D2 = zeros(2, Unknowns);
D2(1,1)     = 3*X(1)^2;
D2(1,2)     = 2*X(1);
D2(1,3)     = 1;
D2(2,end)   = 1;
D2(2,end-1) = 2*X(N);
D2(2,end-2) = 3*X(N)^2;
% Stack the equations into a matrix
H = [H0; K1; K2; D2];
% Solve the system
ABCs = inv(H)*B;
% Plot the splines
for idx = 1:n
    row = 1 + (idx-1)*4;
    a = ABCs(row);
    b = ABCs(row+1);
    c = ABCs(row+2);
    d = ABCs(row+3);
    xspline = linspace(X(idx), X(idx+1), 100);
    yspline = a*xspline.^3 + b*xspline.^2 + c*xspline + d;
    p3 = plot(xspline, yspline, 'b-');
end
 
% PLOT ALL METHODS
xlim([0 1.05*max(X)])
ylim([0 1.05*max(Y)])
legend([p0,p1,p2,p3], 'Measured Data', 'Linear Least Squares', ...
    'Quadratic Least Squares', 'Cubic Splines')
legend('Location', 'southeast')
 
% LOSS OF ENERGY
% Linear Least Squares Method
LLS = @(x) coefL(1) + x*coefL(2);
InLLS = integral(LLS, X(1), X(N));
% Quadratic Least Squares Method
QLS = @(x) coefQ(1) + x*coefQ(2) + x.^2*coefQ(3);
InQLS = integral(QLS, X(1), X(N));
% Cubic Splines Method
InCS = 0;
for i = 1:n
    a = ABCs(4*(i-1) + 1);
    b = ABCs(4*(i-1) + 2);
    c = ABCs(4*(i-1) + 3);
    d = ABCs(4*(i-1) + 4);
    fun = @(x) a*x.^3 + b*x.^2 + c*x + d;
    In = integral(fun, X(i), X(i+1));
    InCS = InCS + In;
end
% Hooke's Law
InH = .5 * X(N) * Y(N);
% Print values
fprintf('The total energy lost calculated from Linear LS is %5d \n', InLLS')
fprintf('The total energy lost calculated from Quadratic LS is %5d \n', InQLS')
fprintf('The total energy lost calculated from Cubic Splines is %5d \n', InCS')
fprintf('The total energy lost calculated from Hooke''s Law is %5d \n', InH')
 
% CALCULATING ERRORS
abserr_LLS = 0;
abserr_QLS = 0;
for i = 1:length(X)
    abserr_LLS = abserr_LLS + abs(Y(i)-LLS(X(i)));
    abserr_QLS = abserr_QLS + abs(Y(i)-QLS(X(i)));
end
abserr_LLS;
abserr_QLS;
relerr_LLS = abserr_LLS / sum(Y);
relerr_QLS = abserr_QLS / sum(Y);
