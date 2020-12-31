function [Prob linModel dOri MotionDir] = modelFunctionNielsenLempel2021(Params,Data,V1Resp)
% Returns the -log(likelihood) of empirical data given by Data from a model
% implementation with V1 responses given by V1Resp and parameters given by
% Params. 

% Input variables:

% Params is a vector containing the values for paramters: [I Ke Ki Tpss MxResp]
% Data is a structure containing the following field: 
% BResp: the baseline firing rate of the neuron calculated as the 
% spikes/second fiered during presentation of blank stimulus.
% Data: Vector containing the total number of spikes fired by the neuron
% across all presentations for each stimuli. Data stored in psition X of this vector
% must correspond to spikes elicited by a stimulus of dOri and motion
% direction given in the X position of the output variables dOri and
% MotionDir respectively.
% Reps: Vector containing the number of times that a each stimulus
% condition was presented to the animal. Same stimuli order as Data field.
% V1Resp: V1 responses to be used for the model. For this input variable
% responses must be organized in a 2D matrix where where V1Resp(x,y)
% corresponds to the response of the model V1 neuron to a stimuli composed
% of 2 gratings with motion directions 22.5*x and 22.5*y deg. Positions along
% the diagonal indicate responses to gratings of direction 22.5*x deg. 

% Output variables:

% Prob: sum of -log(likelihood) of empirical data given by Data assuming a
% poisson distribution of spikes for each condition with event frequency
% equal to the response of the model PSS neurons.
% linmodel: responses for the modeled PSS neuron for all stimuli
% conditions. Same stimuli order as Data.Data.


Inh=Params(1);
Width=Params(2);
WidthI=Params(3);
MtNL=Params(4);
MXResp=Params(5);


V1Resp=V1Resp./max(V1Resp(:));

%Compute V1 population responses

for N = 1:16
    popV1Resp(:,:,N) = circshift(V1Resp,[N-1 N-1]);
end

% Set MT Weights

MTE = circ_vmpdf([0:1:15]*(2*pi/16),14*pi/16,2^Width);
MTE = MTE./max(MTE);

MTI = circshift(circ_vmpdf([0:1:15]*(2*pi/16),14*pi/16,2^WidthI),8);
MTI = MTI./max(MTI);

MTW = MTE - MTI*Inh;

% Get MT Resp

MtModel = zeros(16,16);
 
for i = 1:16
    MtModel = MtModel + popV1Resp(:,:,i).*MTW(i);
end

if max(MtModel(:))>0
MtModel=MtModel./(max(MtModel(:)));
else
MtModel = zeros(16,16);    
end
MtModel(MtModel<0)=0;
MtModel = MtModel.^MtNL;
% MtModel = MtModel - MtNL;

if max(MtModel(:))>0
MtModel=MtModel./(max(MtModel(:)));
else
MtModel = zeros(16,16);    
end

% Allign to pref ori

for ori = 1:16
TCunit(ori) = MtModel(ori,ori);
end

mx = find(TCunit == max(TCunit));

MtModel = MtModel([mx:end 1:mx-1],[mx:end 1:mx-1]);
linModel=[];
for i = 1:16
    linModel=[linModel MtModel(i,1:end-(i-1))];
end
linModel([9 25 40 54 67 79 90 100])=[];

linModel=linModel*MXResp;
linModel=linModel+Data.BResp;
% Implement minimum baseline response frequency to prevent expected
% response form being 0 which would result in log(likelihood) = -inf.
linModel(linModel<.01)=.01;

% For every empirical response calculate the probability of encoutered
% events (spikes) given a poisson distribution with event frequency equal
% to model response
for i =1:length(linModel)
    Probs(i)=poisspdf(round(Data.Data(i)),(linModel(i)*Data.Reps(i)*.15));
end
% Sum the -log(likelihood) of all empirical responses
Prob = -sum(log(Probs));

% Make dOri and MotionDir output variables
dOri=[0:22.5:157.5 157.5:-22.5:22.5 0:22.5:157.5 157.5:-22.5:45 0:22.5:157.5 157.5:-22.5:67.5 0:22.5:157.5 157.5:-22.5:90 0:22.5:157.5 157.5:-22.5:102.5 0:22.5:157.5 157.5 135 0:22.5:157.5 157.5 0:22.5:157.5 0:22.5:157.5 0:22.5:135 0:22.5:112.5 0:22.5:90 0:22.5:67.5 0:22.5:45 0 22.5 0];
MotionDir= [0:11.25:78.75 281.25:11.25:348.75 [0:11.25:78.75 281.25:11.25:337.5]+22.5 [0:11.25:78.75 281.25:11.25:326.25]+45 [0:11.25:78.75 281.25:11.25:315]+67.5 [0:11.25:78.75 281.25:11.25:303.75]+90 [0:11.25:78.75 281.25:11.25:292.5]+112.5 [0:11.25:78.75 281.25]+135 [0:11.25:78.75]+157.5 [0:11.25:78.75]+180 [0:11.25:67.5]+202.5 [0:11.25:56.25]+225 [0:11.25:45]+247.5 [0:11.25:33.75]+270 [0 11.25 22.5]+292.5 [0 11.25]+315 337.5];
MotionDir(MotionDir>=360)=MotionDir(MotionDir>=360)-360;
