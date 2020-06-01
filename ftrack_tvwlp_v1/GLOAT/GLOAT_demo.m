function [] = GLOAT_demo()

[wave,Fs]=wavread('arctic_awb_a0001.wav');

F0min=80;
F0max=240;

%% Pitch tracking using a method based on the Summation of the Residual Harmonics
[f0,VUVDecisions,SRHVal] = SRH_PitchTracking(wave,Fs,F0min,F0max);

VUVDecisions2=zeros(1,length(wave));
HopSize=round(10/1000*Fs);
for k=1:length(VUVDecisions)
    VUVDecisions2((k-1)*HopSize+1:k*HopSize)=VUVDecisions(k);    
end

figure
subplot(211)
plot(f0)
xlabel('Frame Index')
ylabel('F0 (Hz)')
title('Pitch Tracking with the SRH algorithm')
subplot(212)
plot(SRHVal)
hold on
plot(0.07*VUVDecisions,'r')
xlabel('Frame Index')
ylabel('SHR (F0)')
hold off
pause(0.00001)

%% Estimation of the mean pitch value
f0_tmp=f0.*VUVDecisions;
pos= f0_tmp~=0;
f0_tmp=f0_tmp(pos);
F0mean=mean(f0_tmp);


%% Oscillating Moment-based Polarity Detection
[Polarity] = OMPD_PolarityDetection(wave,Fs,F0mean,VUVDecisions);

%% Speech Event Detection using the Residual Excitation And a Mean-based Signal
[gci,MeanBasedSignal] = SEDREAMS_GCIDetection(Polarity*wave,Fs,F0mean);

% Optimized version of the code in terms of computation time
%[gci,MeanBasedSignal] = SEDREAMS_GCIDetection_ComputationalPerformanceOptimized(Polarity*wave,Fs,F0mean);

res = GetLPCresidual(wave,round(25/1000*Fs),round(5/1000*Fs),round(Fs/1000)+2);
res=Polarity*res/max(abs(res));

figure
plot(res)
hold on
plot(gci,res(gci),'xr')
plot(MeanBasedSignal,'m')
plot(VUVDecisions2,'k')
legend('Residual excitation','GCI locations','Mean-based signal','V-UV decisions')
title('GCI detection with the SEDREAMS algorithm')
pause(0.00001)

%% Glottal Source estimation using the Complex Cepstrum
[GlottalSource] = CCD_GlottalFlowEstimation(Polarity*wave,Fs,gci,f0,VUVDecisions);

figure
subplot(211)
plot(wave/max(abs(wave)))
hold on
plot(VUVDecisions2,'r')
xlabel('Time (samples)')
ylabel('Speech Signal')
title('Glottal Source Estimation with the Complex Cepstrum-based Decomposition')
subplot(212)
plot(GlottalSource)
xlabel('Time (samples)')
ylabel('Glottal Source Estimate')
pause(0.00001)