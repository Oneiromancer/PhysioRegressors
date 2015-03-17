%directory of acqknowledge txt file
Raw_Dir = 'D:\PHLEM & RVHR\PhysioRegressors\TestData';
Raw_PhysioTxt = {'RestingState_7min_500hz'}; %don't include txt extension
    
TR = 2;

%%
%makes and saves .mat file with dataOut
cutPhysioData(Raw_PhysioTxt, TR)

%load dataout
cutPhysioStruct = ls('*TRs.mat');
load(cutPhysioStruct)
%%

%optional -- load dataOut with PhysioScoreSMG


%%
%Cut off extra time beyond number of TRs 
PhLEM = BIOPAC_to_PhLEM(dataOut)
%%
%make regressors
%once for HR, once for RESP

%HeartRate
hr_field = 'PPOHIGH';
[HRTS names] = get_rate_timeseries(hr_field,PhLEM);

%RespRate
resp_field = 'RESP';
[RespTS names] = get_rate_timeseries(resp_field, PhLEM);

save('HRTS', 'HRTS')
display('Saved Heart Rate Time Series')
save('RespTS', 'RespTS')
display('Saved Respiration Rate Time Series')
    

%%
%Chang RVHRCor

TR = 2;

ana_type = 'rvhr';
plotflag = 0;

% Check for options to reset DEFAULTS
% for v=1:2:length(varargin),
%     eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
% end

% Get the target field of the PhLEM Structure
%dtmp = eval(sprintf('PhLEM.%s',param));

% Setup TR markers
tr_sampling = round(PhLEM.Time(end)/TR);
tr_timestamps = [0:(TR):PhLEM.Time(end)-(TR/2)]/mean(diff(PhLEM.Time));

hr_fields = {'ECG','PPOHIGH'};
rv_fields = {'RESP','PPOLOW'};
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names = {'RVRRF'};
%compute RV regressor 
D = PhLEM.RESP;
Time = PhLEM.Time;
    
num_tr = round(Time(end)/TR);
tr_vec = floor(Time/TR)+1;

% compute RV vector according to Chang et al:
for i=1:num_tr
  rv(1,i) = std(D.smoothed_data(ismember(tr_vec,[i-1,i,i+1])));
end

% normalize using z-score:
rv = zscore(rv);

% ...convolve with RRF (from Reference 2):
t = 0:TR:28;
RRF = 0.6*t.^2.1.*exp(-t/1.6)-0.0023*t.^3.54.*exp(-t/4.25);
Col = conv(rv,RRF);

% ...trim excess:
Col = Col(1:num_tr)';

RVRRF = struct('names', names, 'RVRRF', Col);
save('RVRRF', 'RVRRF')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute HR regressor 
names = {'HRCRF'};
%compute RV regressor 
D = PhLEM.PPOHIGH;
Time = PhLEM.Time;


num_tr = round(Time(end)/TR);
beats = Time(D.events==1);
tr_vec = floor(beats/TR)+1;

% compute HR vector according to Chang et al:
for i=1:210  %num_tr changed to exact number due to it thinking more than 210 TRs, likekly because end of PPO is not cutoff
  hr(1,i) = 60/mean(diff(beats(ismember(tr_vec,[i-1,i,i+1]))));
end

% normalize using z-score:
hr = zscore(hr); %hr-mean(hr);

% ...convolve with CRF (from Reference 1):
t = 0:TR:28;
CRF = 0.6*t.^2.7.*exp(-t/1.6) - 16/sqrt(2*pi*9)*exp(-0.5*(t-12).^2/9);

% And in case I need the CRF's time derivative:
CRF_TD = (1.62*t.^1/7 - (3/8)*t.^2.7).*exp(-t/1.6) + ...
	 (t./9 + 8/3).*exp(-(t-12).^2/18);

Col = conv(hr,CRF);

Col2 = conv(hr,CRF_TD);


% ...trim excess:
Col = Col(1:num_tr)'; 
Col2 = Col2(1:num_tr)';

HRCRF = struct('names', names, 'HeartRateConv', Col, 'HeartRateConvwithTimeDerivative', Col2);
save('HRCRF', 'HRCRF')