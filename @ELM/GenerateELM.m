function elm=GenerateELM(R,N,activeFcn,S)
%Train ELM with input dimension as R*Q, N hidden neurons with acitivation function as activeFcn, output dimension as S
    % Randomly Generate the Input Weight Matrix
    elm.weight{1} = rand(N,R) * 2 - 1;
    elm.weight{2} = rand(S,N) * 2;
    % Randomly Generate the Bias Matrix
    elm.bias{1} = rand(N,1);
    elm.bias{2} = zeros(S,1);
    elm.layerNum = 2;
    elm.activeFcn = {activeFcn  'purelin'};
    elm.trainingError = [];
    elm=ELM(elm.weight,elm.bias,elm.activeFcn);
end
