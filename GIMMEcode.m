%Integration of Gene Expression Data into Lachancea kluyveri model

initCobraToolbox()

cd D:/ScientificReports/

model=readCbModel('iPN730_etac_rich')


cd D:/ScientificReports/GeneExpressionIntegration/
expr=readtable('NormalizedAbsoluteExpression.csv')

expressionData.gene=table2array(expr(:,5))

exp=table2array(expr(:,4))
exp_val=str2double(exp)
expressionData.value=exp_val

[expressionRxns, parsedGPR] = mapExpressionToReactions(model, expressionData)

p = 0:0.25:1;
y1 = quantile(expressionRxns,p);
z_i = [p;y1]

subplot(1,2,1)
hist(expressionRxns,50)
subplot(1,2,2)
boxplot(expressionRxns)

thr=-5:0.4:5

for i=1:length(thr)
    tissueModel = GIMME(model, expressionRxns, thr(i))
    rxn_no(i)=length(tissueModel.rxns)
    met_no(i)=length(tissueModel.mets)
    sol=optimizeCbModel(tissueModel)
    gr(i)=sol.f
    etac_bool(i)=length(tissueModel.rxns(ismember(tissueModel.rxns,'AAT')))
    etoh_bool(i)=length(tissueModel.rxns(ismember(tissueModel.rxns,'R_ALCD2ir_c0')))
end

plot(thr,rxn_no,thr,met_no)
hold on
bar(etac_bool)

tissueModel = GIMME(model, expressionRxns, -0.4)
etac_bool=length(tissueModel.rxns(ismember(tissueModel.rxns,'AAT')))
etoh_bool=length(tissueModel.rxns(ismember(tissueModel.rxns,'R_ALCD2ir_c0')))

tissueModel=changeRxnBounds(tissueModel,'EX_cpd00027_e0',-2.28,'l')

%Aerobic conditions
tissueModel=changeRxnBounds(tissueModel,'EX_cpd00007_e0',-10,'l')
sol=optimizeCbModel(tissueModel)
printFluxVector(tissueModel,sol.x,'true','true')


%Semoaerobic conditions
tissueModel=changeRxnBounds(tissueModel,'EX_cpd00007_e0',-3,'l')
sol=optimizeCbModel(tissueModel)
printFluxVector(tissueModel,sol.x,'true','true')

%Anaerobic Conditions
tissueModel=changeRxnBounds(tissueModel,'EX_cpd00007_e0',-0.25,'l')
sol=optimizeCbModel(tissueModel)
printFluxVector(tissueModel,sol.x,'true','true')

%Check Ethanol production
etoh=sol.x(ismember(model.rxn,'EX_cpd00363_e0'))
if (etoh>0)
    fprintf("Ethanol is produced:")
end

%Check ethyl acetate production
etac=sol.x(ismember(model.rxn,'EX_cpd00363_e0'))
if (etoh>0)
    fprintf("Ethanol is produced: %d")
end


e=expressionRxns
sol_norm=(sol.x-mean(sol.x))/std(sol.x)
exp_norm=(e-nanmean(e))/nanstd(e)
hist(sol_norm)
hist(exp_norm)
