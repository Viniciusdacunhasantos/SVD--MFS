function Y = sphemfr(l_max, m_max, th,lam)
%Adapted from MYSPHHARM   Main subroutine for computing the spherical harmonics.
%
% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Below is the implementation the Modified Forward Row (MFR) method 
% described in the Holmes and Featherstones paper (2002)
% It computes P^m_n/u^m by a stable recurrence to avoid numerical errors
% near the poles

abs_m_max = abs( m_max );

% Make lam row
%lam = lam.';
% Make th vector
%th = th.';

% Precompute vector cos(th)
p = length(th);
CosTh = cos(th);

%% Compute P_m_max^m_max / u^m_max

% Initialize P^0_0 / u^0
Pold = ones(p, 1);

% Compute P^m_m/u^m
for m = 1:l_max
    % Compute Pm^m with Pm-1^m-1
    Pold = sqrt( (2*m+1)/(2*m-(m==1)) ) * Pold;
end

% Initialize the recurrence (Pm^m-1 does not exist)
Poldold = zeros(p, 1);

%% Compute P^m_l / u^m with the recurrence formula, m_max+1 <= l <= l_max
for m = l_max-[1:(l_max-abs_m_max)]
    gnm =2*(m+1)/sqrt((l_max-m)*(l_max+m+1));
    hnm=sqrt((l_max+m+2)*(l_max-m-1)/((l_max-m)*(l_max+m+1)));
    % Compute the normalized associated legendre polynomial P^m_l/u^m
    Pl = 1/sqrt(2^(m==0))*(gnm*CosTh.*Pold - hnm*((sin(th)).^2).*Poldold); %(1/sqrt(2^(m==0)))*
    %Pl = (gnm.*CosTh.*Pold);
    
    % Update the polynomials for the recurrence
    Poldold = Pold;
    Pold = Pl;

end

% Normalize the polynomial and recover associated Legendre polynomials
Pold = ((-1)^m_max)*(sin(th).^abs_m_max).*Pold/sqrt(4*pi); %(-1)^abs_m_max*
Y=Pold.*exp(1i*m_max*lam)/(sqrt(2^(m_max~=0)));
end
    