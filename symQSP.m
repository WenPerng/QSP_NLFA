%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% symQSP.m
%--------------------------------------------------------------------------
% Input:
%   x   -- coordinates from -1 to 1, e.g. x = linspace(-1,1,101).
%   phi -- phase factors (phi_0,...,phi_d).
%
% Output:
%   curve -- a function f(x) from x=-1 to 1, where
%               f(x) = Im{U11},
%               U11  = [U(x)]_{11},
%               U(x) = exp(i(phi_d)Z) * V(x) * exp(i(phi_{d-1})Z) * V(x) ...
%                      * exp(i(phi_1)Z * V(x) * exp(i(phi_0)Z) * V(x) ...
%                      * exp(i(phi_1)Z * ... * V(x) * exp(i(phi_d)Z),
%               V(x) = [x, i sqrt(1-x^2); i sqrt(1-x^2) , x]
%                    = exp(i acos(x) X),
%            which is simply the symmetric QSP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function curve = symQSP(phi,x)
    W = @(t) [t,1i*sqrt(1-t^2);1i*sqrt(1-t^2),t];
    U = @(t) QSP_mult(W(t),phi);
    curve = arrayfun(@(t) [1,0] * U(t) * [1;0],x);
    curve = imag(curve);
end

function U = QSP_mult(W,phi)
    Z = [1,0;0,-1];
    U = expm(1i * phi(1) * Z);
    for i = 2:length(phi)
        U = expm(1i * phi(i) * Z) * W * U * W * expm(1i * phi(i) * Z);
    end
end