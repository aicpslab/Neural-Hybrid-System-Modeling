function [partitions,ELMs] = MergePatitions(obj,ELMs1,gamma,mu,coeff,inputdimension)

xs=obj.input;
t=obj.output;
segmentIndex=obj.intervals;
inputspace=segmentIndex;
lengthSegment=size(obj.intervals,2);
if nargin <4
for k = 1:1:size(obj.intervals,2) 
    [TrainDatainput{k},TrainDataoutput{k}]=obj.Dataselect(xs,t,obj.intervals{k},2);                      
end
elseif nargin==5
 dimension = size(obj.intervals{1},1) ;   
for k = 1:1:size(obj.intervals,2) 
    [TrainDatainput{k},TrainDataoutput{k}]=obj.Dataselect(xs,t,obj.intervals{k},dimension,mu,coeff);                      
end
elseif nargin==6
 dimension = size(obj.intervals{1},1) ;   
    for k = 1:1:size(obj.intervals,2) 
        [TrainDatainput{k},TrainDataoutput{k}]=obj.Dataselect(xs,t,obj.intervals{k},dimension,mu,coeff,inputdimension);                      
    end    
end
inputDataSegmented=TrainDatainput;
outputDataSegmented=TrainDataoutput;
fprintf('Starting Merging redundant partitions, there are %d',k);
fprintf('partitions')
m1 = 1;
k = 0;

while m1 <  lengthSegment
    m2 = m1 + 1;
    while m2 <= lengthSegment
        inputData = [inputDataSegmented{m1}{1},inputDataSegmented{m2}{1}];
        outputData = [outputDataSegmented{m1}{1},outputDataSegmented{m2}{1}];
        %ELMNetwork = trainELMLipridge(ELMs1,inputData,outputData);
        ELMNetwork = trainELM(ELMs1,inputData,outputData);
        k = k + 1;
        if ELMNetwork.trainingError < gamma
            inputDataSegmented{m1}{1} = inputData;
            inputDataSegmented(m2) = [];
            outputDataSegmented{m1}{1} = outputData;
            outputDataSegmented(m2) = [];
            segmentIndex{m1} = [segmentIndex{1,m1};segmentIndex{1,m2}];
            segmentIndex(:,m2) = [];
            lengthSegment = lengthSegment - 1;
        else
            m2 = m2 + 1;
        end
    end
    m1 = m1 + 1;
end

clc
for i = 1:size(inputDataSegmented,2)
    if(~isempty(TrainDatainput{1,i}))
         ELMs(i)=trainELMLipridge(ELMs1,inputDataSegmented{i}{1},outputDataSegmented{i}{1});
       %  ELMs(i)=trainELM(ELMs1,inputDataSegmented{i}{1},outputDataSegmented{i}{1}); 
       inputspace1{1,k}=inputspace{1,i};
    end
end
obj.intervals=[];
obj.intervals=segmentIndex;
partitions = obj;
end