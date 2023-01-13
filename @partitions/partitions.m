classdef partitions
    % partitions class
    %
    % Syntax:
    %    object constructor: Obj = ffnetwork(w,b,a)
    %    copy constructor: Obj = otherObj
    %
    % Inputs:
    %    input1 - weight, cell{matrix}
    %    input2 - bias, cell{vector}
    %    input3 - activation function, cell{string}
    %
    % Outputs:
    %    Obj - Generated Object
    %
    
    % Author:       Yejiang Yang
    % Written:      02/25/2019
    % Last update:  09/13/2019
    
%------------- BEGIN CODE --------------
    
    properties
        intervals = {};
        input = {};
        output = {};
    end
    
    methods
        %% class constructor
        function obj = partitions(intervals,input,output)
            obj.intervals = intervals;
            obj.input = input;
            obj.output=output;
        end
        
        %% methods in seperate files
        
        % generate ELM model with input matrix P and output target matrix T
        y = ME(obj,length,epsilon,dimension,mu,coeff)
        
        % Merge redundant partitions with tolerance as gamma
        [partitions,ELMs] = MergePatitions(obj,ELM1,gamma,mu,coeff,inputdimension)

        % Spliting computation for each intersection 
  

        % Combine the output reachable set for one time step
        
 
%         % generate numInput random outputs for an input interval inputInterval
%         y = output(obj,input)


%         % generate output set with nnv as reachability computational method 
%       
%         y = SetNNV(obj,inputset)
%         
%         % generate output set with interval as reachability computational method 
%         y = SetInterval(obj,inputset)

        %

        % pre-generate equations for computing output interval
%         generateFcn(obj)
%         
%         % compute the output set for a single interval input
%         outputIntvl = outputSetSingle(obj,inputIntvl)
%         
%         % compute the output set for an interval using UNIFORM PARTITION method
%         [outputIntvl_array,inputIntvl_array] = outputSetUniform(obj,inputIntvl,tol) 
%         
%         % compute the output set for an interval using SPECIFICATION-GUIDED method method
%         [outputIntvl_array,inputIntvl_array] = outputSetSpec(obj,inputIntvl,unsafeIntvl,tol) 
%         
%         % compute the output set for an interval using SIMULATION-GUIDED method
%         [outputIntvl_array,inputIntvl_array] = outputSetSim(obj,inputIntvl,options)
%         
%         % compute the input set for the target interval set using BACKWARD COMPUTATION method
%         [outputIntvl_array,inputIntvl_array] = outputSetBack(obj,inputIntvl,targetIntvl,tol) 
%         
%         % compute the input set for a feedforward network
%         [outputIntvl_array,inputIntvl_array] = outputSet(varargin)
% 
%         % compute the output set for NARMA model using UNIFORM PARTITION
%         % method
%         outputIntvl_array = outputSetNARMAUniform(obj,initialIntvl,inputIntvl,step,tol)
%         
%         % compute the output set for NARMA model using SIMULATION-GUIDED
%         % method
%         outputIntvl_array = outputSetNARMASim(obj,initialIntvl,inputIntvl,step,options)
%         
%         % compute the output set for NARMA model using SIMULATION-GUIDED
%         % method
%         outputIntvl_array = outputSetNARMA(varargin)
        
        
    end
    methods(Static)

% Selecte Data within one partitions
[input,output]=Dataselect(input,output,Xtemp2,dimension,mu,coeff,inputdimension)
%
flag=ifin(data,space,dimension)

% compute the distance between two neural networks
delta_entropy = cal_entropy(Xtemp1,Xtemp2,xs,t,dimension)
   
error= MeanSquare(Input)

[splitingset,t] = SplitingReach(intervals,ELM,setinput,Reachmethod,ubound)

InterSet=SetIntersect(setinput,partitions)

intervalplot(I,fillType,color)

y= CombineReach(setinput)
           
%         disty = distOutputSingle(ffnn1,ffnn2,x)
%         
%        % compute the distance of random outputs of two nerual networks for an input interval inputInterval
%         disty = distOutputRandom(ffnn1,ffnn2,inputIntvl,num)
%         
%         % compute the distance between two neural networks for one single
%         % input interval
%         distOutputIntvl = distOutputSetSingle(ffnn1,ffnn2,inputIntvl)
%         
%         % compute the distance between two neural networks for a single
%         % input interval using UNIFORM partition method
%         [distOutputIntvl_array,inputIntvl_array] = distOutputSetUniform(ffnn1,ffnn2,inputIntvl,tol)
%         
%         % compute the distance between two neural networks for a single
%         % input interval using simulation partition method
%         [outputIntvl_array,inputIntvl_array] = distOutputSetSim(ffnn1,ffnn2,inputIntvl,options)
%         
%         %compute the distance between two neural networks for a single
%         % input interval 
%         [outputIntvl_array,inputIntvl_array] = distOutputSet(varargin)
%         
%         % compute the distance between two NARMA model using uniform
%         % partition method
%         outputIntvl_array = distOutputSetNARMAUniform(ffnn1,ffnn2,initialIntvl,inputIntvl,step,tol)
%         
%         % compute the distance between two NARMA model using uniform
%         % simulation-guided method
%         outputIntvl_array = distOutputSetNARMASim(ffnn1,ffnn2,initialIntvl,inputIntvl,step,options)
%         
%         % compute the distance between two NARMA model
%         outputIntvl_array = distOutputSetNARMA(varargin) 
     end
end
%------------- END OF CODE --------------
    
    
    
    
