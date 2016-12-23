clear;

data = randi([0 1],96,1);
hModulator = comm.QPSKModulator('BitInput',true);
hModulator.PhaseOffset = 0;
modData = step(hModulator, data);
% scatterplot(modData)