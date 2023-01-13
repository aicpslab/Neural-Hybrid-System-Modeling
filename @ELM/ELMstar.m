function starout=ELMstar(star,ELM)
[N,S]=size(ELM.weight{1,1});
R=size(ELM.weight{1,2},1);
inputweight=ELM.weight{1,1};
outputweight=ELM.weight{1,2};
inputbias=ELM.bias{1,1};
outputbias=ELM.bias{1,2};
L1=LayerS(inputweight,inputbias,'poslin');
L2=LayerS(outputweight,outputbias,'purelin');
ELMnnv=FFNNS([L1 L2]);
cores=4;
[starout,~]= ELMnnv.reach(star,'exact-star',cores);
end