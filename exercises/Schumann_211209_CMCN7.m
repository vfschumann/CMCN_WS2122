% Work of Vic-Fabienne Schumann and Stefan Bauermeister 7816664

%ROOM 

% Task1
% emodel=readCbModel
% emodel=creategrRulesField(emodel)
% Task2
changeCobraSolver('glpk')
w=optimizeCbModel(emodel)

%Task3

wu=(1:length(emodel.rxns))
wl=(1:length(emodel.rxns))
wu=wu'
wl=wl'

for i = 1:(length(emodel.rxns))
    wu(i)=w.x(i)+0.03*abs(w.x(i))+0.001
    wl(i)=w.x(i)-0.03*abs(w.x(i))-0.001
end

% Task4
% for max constraint
I = zeros([length(emodel.rxns),length(emodel.rxns)])
        for x = 1:length(emodel.rxns)
            for y = 1:length(emodel.rxns)
                if x==y
                I(x,y)=1;
                end
            end
        
        end



Au1 = I
Au2 = I
for x = 1:length(emodel.rxns)
            for y = 1:length(emodel.rxns)
                if x==y
                Au2(x,y)=-emodel.ub(x)+wu(x);
                end
            end
end

Au= [Au1 Au2]

% for min constraint

Al1 = -I
Al2 = I
for x = 1:length(emodel.rxns)
            for y = 1:length(emodel.rxns)
                if x==y
                Al2(x,y)=emodel.lb(x)-wl(x);
                end
            end
end

Al=[Al1 Al2]

bu=wu
bl=-wl

% f vector

f=(1:2*length(emodel.rxns))

for x = 1:2*length(emodel.rxns)
    if x <= length(emodel.rxns)
        f(x)=0
   
    else
        f(x)=1
    end
end

% equality Constraints
 o =zeros([length(emodel.mets),length(emodel.rxns)])

Aeq=[emodel.S o]
beq=(1:length(emodel.mets))
beq=beq'
beq(1:length(emodel.mets))=0


yub = (1:length(emodel.rxns))
ylb =(1:length(emodel.rxns))
yub(1:length(emodel.rxns))= 1
ylb(1:length(emodel.rxns))= 0
yub=yub'
ylb=ylb'

lb=[emodel.lb;ylb]
ub=[emodel.ub;yub]

findRxnsFromGenes(emodel,emodel.genes(32))
% reaction 34 is from gene b2779
lb(34)=0
ub(34)=0

A=[Au;Al]
b=[bu;bl]

intcon = length(emodel.rxns)+1:2*length(emodel.rxns)

[xu,fvalu] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub)

% Task 6 calculate euclideandistance
distance1=(1:length(emodel.rxns))
distance1=distance1'


for i = 1:length(emodel.rxns)
    distance1(i)= (w.x(i)-xu(i))^2
end

euclideandistance1 = (sum(distance1))^(1/2)

%MOMA
emodeldel=emodel
emodeldel.ub(34)=0
emodel.lb(34)=0



[solution1, solution2, totalFluxDiff] = optimizeTwoCbModels(emodel, emodeldel)



distance2=(1:length(emodel.rxns))
distance2=distance2'


for i = 1:length(emodel.rxns)
    distance2(i)= (solution1.x(i)-solution2.x(i))^2
end
euclideandistance2=(sum(distance2))^(1/2)



