% Calculates equivalent linear optical loss and amplitude coupling
% coefficient of a waveguide-coupled optical resonator, based on the given
% cavity Q-factor and extinction ratio
% 
% This function only applies to a resonator coupled to a single waveguide,
% i.e. a resonator without the drop port; and this function only applies to
% resonators that have a finesse >> 1
%
% Author: JJ Hu @ MIT
%
% [loss, k] = Resonator_Loss(Q, extinction, FSR, L, lamda, jj)
% loss:     equivalent linear optical loss (waveguide loss) of resonator in
%           dB/cm
% k:        the amplitude coupling coefficient (the phase term neglected)
% Q:        cavity quality factor
% extinction:   extinction ratio in dB (please use a POSITIVE number!)
% FSR:      free spectral range in nm
% L:        physical length/circumference of resonator in um
% lamda:    approximate resonant wavelength in um (e.g. 1.55)
% jj is a variable that indicates the coupling regime; jj = 0 corresponds
% to undercoupling and jj <> 0 corresponds to overcoupling
%
% Reference: A. Yariv, "Universal relations for coupling of optical power
% between microresonators and dielectric waveguides"

function [loss, k] = Resonator_Loss(Q, extinction, FSR, L, lamda, jj)

ng = lamda.^2./L./FSR*1000;  % ng: group index of resonant mode
A = pi.*ng.*L./Q./lamda;
Tmin = 10.^(-extinction/10);
B = A*sqrt(Tmin);

if jj == 0
    alpha = 1/2*(-B+sqrt(B.^2-4*(A-1)));
    t = (1-A)./alpha;
    loss = -2./L.*log(alpha)*10000*4.3429;
    k = sqrt(1-t.^2);
else
    t = 1/2*(-B+sqrt(B.^2-4*(A-1)));
    alpha = (1-A)./t;
    loss = -2./L.*log(alpha)*10000*4.3429;
    k = sqrt(1-t.^2);
end