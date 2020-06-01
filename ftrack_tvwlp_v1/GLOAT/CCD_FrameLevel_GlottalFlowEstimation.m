function [minPhase,maxPhase] = CCD_FrameLevel_GlottalFlowEstimation(inSignal)

% USAGE:
% [minPhase,maxPhase] = ComplexCepstrumDecomposition(inSignal)
% 
% where inSignal is a frame of voiced speech properly windowed,
% minPhase and maxPhase are respectively the minimum- and maximum-phase
% components of the speech signal. The mixed-phase separation is achieved
% using the Complex Cepstrum-based algorithm.
% 
% Windowing is critical and should respect the conditions
% detailed in:
% 
% T.Drugman, B.Bozkurt, T.Dutoit, "Causal-anticausal Decomposition of 
% Speech using Complex Cepstrum for Glottal Source Estimation",
% Speech Communication, Volume 53, Issue 6, July 2011, Pages 855-866, 2011
% 
% and
% 
% Thomas Drugman, Baris Bozkurt, Thierry Dutoit,
% "Complex Cepstrum-based Decomposition of Speech for Glottal Source Estimation",
% Proc. Interspeech09, Brighton, U.K, 2009.
% 
% Please refer to these works in your publication if you use this code.
% 
% Code written by Thomas Drugman in TCTS Lab, University of Mons, Belgium.
% 
% Copyright (C) 2000-2011 Thomas Drugman - TCTS Lab
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Matlab version used for developing this code: 7.0
%
% http://tcts.fpms.ac.be/~drugman/      thomas.drugman@umons.ac.be


halfLength=round(length(inSignal)/2);

[xhat]=ComputeComplexCepstrum(inSignal);   

antiCausalPortionCepstrum=xhat;
antiCausalPortionCepstrum(2:round(length(xhat)/2))=0;
antiCausalPortionCepstrum(1)=antiCausalPortionCepstrum(1)/2;
           
causalPortionCepstrum=xhat;
causalPortionCepstrum(round(length(xhat)/2)+1:end)=0;
causalPortionCepstrum(1)=causalPortionCepstrum(1)/2;

minPhase=inverse_CC(causalPortionCepstrum);
maxPhase=inverse_CC(antiCausalPortionCepstrum);
    
maxPhase=diff(maxPhase);
minPhase=filter([1 0],[1 -1],minPhase);

maxPhase=maxPhase(end-halfLength+1:end);
if(max(maxPhase)>max(-maxPhase))
    maxPhase=-maxPhase;    
end

minPhase=minPhase(1:length(inSignal)-length(maxPhase));
if(max(minPhase)>max(-minPhase))
    minPhase=-minPhase;    
end

maxPhase=maxPhase/(max(abs(maxPhase)));
minPhase=minPhase/(max(abs(minPhase)));



function [xhat] = ComputeComplexCepstrum(x)

n=2048;

h = fft(x,n);
[ah] = rcunwrap(angle(h));

logh = log(abs(h))+1i*ah;
xhat = real(ifft(logh));

function [y,nd] = rcunwrap(x)
%RCUNWRAP Phase unwrap utility used by CCEPS.
%   RCUNWRAP(X) unwraps the phase and removes phase corresponding
%   to integer lag.  See also: UNWRAP, CCEPS.

%   Author(s): L. Shure, 1988
%   	   L. Shure and help from PL, 3-30-92, revised

n = length(x);
y = unwrap(x);

nh = fix((n+1)/2);
y(:) = y(:)' - y(nh+1)*(0:(n-1))/nh;


function x = inverse_CC(xhat)

logh = fft(xhat);
h = exp(real(logh)+1i*imag(logh));
x = real(ifft(h));

