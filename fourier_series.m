function [y_pred,y_pred_,y_pred__] = fourier_series(a_,x)
%FOURIER_SERIES.M This function evaluates a Fourier series given the
%Fourier series coefficients a_ and the parameterized value x.

% Input variables
%x - Parameterization variable. In practice, this is the arc length at a
%point of interest. 

% a_ - Fourier series coefficients.

% Output variables
% y_pred - Predicted variable. In practice, this is a fitted or predicted
% value of the x- or y- component of a point of interest. 

% y_pred_ - First deritivate of perdicted variable 

% y_pred__ - Second derivative of perdicted variable 

% concatenates Fourier series coefficients into a*cos and b*sin
a0 = a_(1);
a = a_(2:2:end);
b = a_(3:2:end);

% pre-allocation of variables
y_pred = zeros(length(x),1);
y_pred_ = zeros(length(x),1);
y_pred__ = zeros(length(x),1);

% evaluates y_pred, y_pred_, and y_pred__
for i = 1:length(x)
    for j = 1:length(a)
        if j == 1
            y_pred(i) = a0 + a(j)*cos(j*x(i)) + b(j)*sin(j*x(i));
            y_pred_(i) = -a(j)*sin(j*x(i))*j + b(j)*cos(j*x(i))*j;
            y_pred__(i) = -a(j)*cos(j*x(i))*j*j - b(j)*sin(j*x(i))*j*j;
        else
            y_pred(i) = y_pred(i) + a(j)*cos(j*x(i)) + b(j)*sin(j*x(i));
            y_pred_(i) = y_pred_(i) + -a(j)*sin(j*x(i))*j + b(j)*cos(j*x(i))*j;
            y_pred__(i) = y_pred__(i) + -a(j)*cos(j*x(i))*j*j - b(j)*sin(j*x(i))*j*j;
        end
    end

end

end