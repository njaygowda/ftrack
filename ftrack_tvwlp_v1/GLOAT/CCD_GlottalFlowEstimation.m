function [GlottalSource] = CCD_GlottalFlowEstimation(wave,Fs,gci,f0,VUVDecisions)

% USAGE:
% [GlottalSource] = CCD_GlottalFlowEstimation(wave,Fs,gci,f0,VUVDecisions)
% 
% INPUTS:
%     - wave is the speech signal
%     - Fs is the sampling frequency (Hz)
%     - gci are the Glottal Closure Instant locations (in samples)
%     - f0 is the vector of F0 values (with an hopsize of 10ms)
%     - VUVDecisions is the vector of voiced-unvoiced VUVDecisions (with an hopsize of 10ms)
%    
% OUPUT:
%     - GlottalSource is a signal containing the derivative of the
%     glottal flow only for voiced glottal cycles
%     
% For more information, please read:
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


display('Glottal Flow Estimation with the Complex Cepstrum-based Decomposition')
pause(0.0001)

% Have a VUV decision and F0 estimation for each sample
VUVDecisions2=zeros(1,length(wave));
f02=zeros(1,length(wave));
HopSize=round(10/1000*Fs);
for k=1:length(VUVDecisions)
    VUVDecisions2((k-1)*HopSize+1:k*HopSize)=VUVDecisions(k);
    f02((k-1)*HopSize+1:k*HopSize)=f0(k);
end

% Keep only GCIs for voiced segments

gindex = VUVDecisions2(gci)~=0;
gcis = gci(gindex);

% res = GetLPCresidual(wave,25/1000*Fs,5/1000*Fs,round(Fs/1000)+2);


% For each glottal cycle, apply a specific windowing satisfying the
% conditions explained in the paper mentioned above, and compute the
% complex cepstrum-based decomposition
GlottalSource=zeros(1,length(wave));

% Low-pass filtering because after separation, an irrelevant
% high-frequency noise (>4000Hz) may appear
Wp = 4000/(Fs/2); Ws = 4100/(Fs/2);
Rp = 3; Rs = 60;
[n,Wn] = ellipord(Wp,Ws,Rp,Rs);
[b_l,a_l] = ellip(n,Rp,Rs,Wn);

for k=1:length(gcis)
    
    T0=round(Fs/f02(gcis(k)));
    
    Seg=wave(gcis(k)-round(0.9*T0):gcis(k)+round(0.9*T0));
    Seg=Seg.*blackman(length(Seg));
    
    [minPhase,maxPhase] = CCD_FrameLevel_GlottalFlowEstimation(Seg);
    [maxPhase] = FilterAndNormalize(maxPhase,b_l,a_l);
    
    maxPhase=maxPhase*sqrt(sum(Seg.^2)/sum(maxPhase.^2));
    
    GlottalSource(gcis(k)-length(maxPhase)+1:gcis(k))=GlottalSource(gcis(k)-length(maxPhase)+1:gcis(k))+maxPhase';    
end



function [sig] = FilterAndNormalize(sig,b_l,a_l)

sig=filtfilt(b_l,a_l,sig);

msig1 = max(sig);
msig2 = abs(min(sig));

if(msig1>msig2)
    sig=-sig;
    msig = msig1;
else
    msig = msig2;
end

sig=sig/(msig);