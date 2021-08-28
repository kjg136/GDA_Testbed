function VarNe(sizes,AnalysisTimes,AvgStart)

for ii=1:length(sizes)
    RunDA(sizes(ii),AnalysisTimes,AvgStart,0,0,1);
end