function [a] = fourier_series_fit(x,y,n)
%FOURIER_SERIES_FIT.M This function computes the coefficients of a fourier
%series for the x or y spatial components of a closed shape. This is done
%to a user-specified degree n which controls the number of fourier series
%used to fit the shape. The total number of coefficients, a, is 2*n+1

% Input variables
%x - Parameterization variable. In practice, this is the arc length at a
%point of interest. 

%y - Variable being fitted. In pratices, this is the x- or y- component of
%a point of interest. 

% n - User-specified degree for Fourier series fit. 

% Output variables
% a - Fitted Fourier series coefficients. The total number of coefficients
% is equal to 2*n+1.

% creates a system of linear equations of the form f*a = y where f contains
% the fourier series evaluations at parameterized value x, a is the fourier series
% coefficents, and y is the variable being fitted. 
f = zeros(length(x),n*2);
for i = 1:length(x)
    for j = 1:n
        f(i,2*j-1) = cos(j*x(i)); % odd columns are cos
        f(i,2*j) = sin(j*x(i)); % even columns are sin
    end
end

% a column of ones appended to the left side of f to account for the
% variable a0
left = ones(length(x),1); 
f = [left,f];

% solves system of linear equations for a
a = f\y;


end