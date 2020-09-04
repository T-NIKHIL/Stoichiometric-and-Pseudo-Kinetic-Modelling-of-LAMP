%                             _ _ _ 
% |          /\     |\    /| |     |
% |         /  \    | \  / | |     |
% |        /_ _ \   |  \/  | |_ _ _|
% |       /      \  |      | |
% |_ _ _ /        \ |      | |

% THIS IS A PROGRAM THAT ANALYZES THE EXPONENTIAL PHASE OF LAMP 

% NAME OF PROGRAM : LAMP_Analyzer_v1

% CAPABILITIES OF THIS PROGRAM :
% 1. IT GENERATES A NETWORK DISPLAYING ALL THE POSSIBLE PRODUCT AMPLICONS THAT CAN BE FORMED FROM A REACTANT AMPLICON.
% 2. IT GENERATES PLOTS OF :
% 2.1 NUMBER OF SS  AMPLICONS GENERATED VS TIME
% 2.2 NUMBER OF PDS AMPLICONS GENERATED VS TIME
% 2.3 NUMBER OF SL  AMPLICONS GENERATED VS TIME
% 2.4 NUMBER OF T   AMPLICONS GENERATED VS TIME
% 2.5 NUCLEOTIDES USED VS TIME

% ASSUMPTIONS:
% 1 We assume primers and nucleotides do not get depleted as we want to find out how many primers and nucleotides are used.
% 2 numOfCopiesPDS/T = (numOfCopiesReactant - numOfCopiesSL)/n .('n' stands for sum of the number of child PDS and T that are formed)
% 3 The probability to form SL from SS or PDS at any given time in the simulation remains the same as the reaction is the same. 
%   Also the rate constants do not vary as it is isothermal reaction.
% 4 We assume all primer binding sites have equal probability for a primer to come and anneal and hence equal concentrations 
%   of each child PDS and T are formed.
% 5 SS can get formed from a PDS at different times due to the strand displacement taking different amounts of times to complete
%   For ex: SS formed from child PDS strand disp will take a smaller time to from than compared to SS formed from a SL strand disp
%   Here we only consider the strand displacement time for SL as a majority of PDS goes to form SL and displaces SS at the same time
% 6 Finally we assume that the child PDS formed from SS or PDS cannot form the remaining larger length child PDS
% 7 T cannot participate in any further reactions
% 8 The anti sense strand of DNA also produces the same amplicons but at slight time differences. For example if you consider sense strand DNA
%    pathway 1 amplicons are formed first and then pathway 2 form but in the case of anti sense first pathway 2 form and then pathway 1 so there is a
%    slight difference in the time the amplicons form (around 3.9 seconds) which is low and we can assume that th amplicons from both sense and
%    anti sense form at the same time
% 9 The dataMat variable stores the time at which a product is formed. This product can react further to form products at larger times, during this
%   transition phase the copies of the reactant are not populated in the time intervals between its inception and degradation

%NOTE 1 :
% 1) PRIMERS : 3' B2 - B1c ; 3' F2 - F1c
% 2) ALL SEQUENCES HAVE 3' AT THE RIGHT MOST END OF ARRAY

%NOTE 2 :
% 1) PDSa  IS PDS FROM SS
% 2) PDSb  IS PDS FROM SL(a or b)
% 3) PDSc  IS PDS FROM PDS(a or b or c)
% 4) SLa   IS SL  FROM SS
% 5) SLb   IS SL  FROM PDS(a or b or c)
% 6) SSa   IS SS  FROM PDSa
% 7) SSb   IS SS  FROM PDSb
% 8) SSc   IS SS  FROM PDSc
% 9) TPDS1 IS T   FROM PDS AND THE 1 REPRESENTS THE INDEX TO ACCESS IN pdsobjArray
% 10)TSS2  IS T   FROM SS  AND THE 2 REPRESENTS THE INDEX TO ACCESS IN ssobjArray

%VARIABLE NOMENCLATURE
% 1) The structure of SL(a or b) is represented using a matrix called ds.
% 2) The structure of SS is represented using an array called ss.
% 3) The structure of PDS is represented using 3 components
% 3.1) ss - this stores the sequence of PDS which will be the future SS amplicon.
% 2.2) ls - this stores the loop sequence of PDS where there exist potential sites for attachement of primers.
% 3.3) es - this stores the sequence of PDS which will be part of the future SL amplicon.

% Examples :
% 1) SLdsPDS  - this variable stores the ds matrix of SL which is formed from PDS.
% 2) PDSssPDS - this variable stores the ss array of PDS which is formed from PDS.
% 3) PDSlsPDS - this variable stores the ls array of PDS which is formed from PDS.
% 4) PDSesPDS - this variable stores the es array of PDS which is formed from PDS.

%------------------------------------------------------------------------------------------------------------------------

% THINGS TO DO BEFORE RUNNING : 
% 1) ENTER ALL INPUT DATA
% 2) MAKE SURE TO RENAME FILE NAMES BEFORE RUNNING THEM TO PREVENT OVERWRITING

%------------------------------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------------------------------

% USER INPUT

% CAUTION : DO NOT RUN FOR MORE THAN 1 MINUTE !!

%Enter the simulation time for which reaction network is to be generated
simTime         = 25;
%Specify the intervals at which to take the readings
timeInt         = 1;
%Enter rate of nucleotide incorporation by the enzyme
b               = 120;
%Specify the initial number of DNA copies
copiesOfDNA     = 100;
%Specify the probability of forming SL (OBTAINED FROM simpleKineticModel.m)
%When ke is 1e8 slformprob = 0.9941
%When ke is 1e7 slformprob = 0.9433
slformprob      = 0.67;

% DO NOT CHNAGE THE VALUES BELOW !

%Specify number of columns to use to store the amplicons
colNum          = 50000;  
%Specify number of rows to use to store the amplicons
%For any simTime less than 1 minute use 50 rows
rowNum          = 50;  
%Specify storage capacity of object arrays (ex: 100 can store 100 amplicons)
capacity        = 50000; 

%Specify starting amplicon sequence
startingSeq     = ["B2","F2c"];

%Specify the type of species (SS,SL,PDS OR T)
type            = "SS";
%Specify the length of bpatch (B1-B2c-B1c/B1-B2-B1c)
B               = 95;
%Specify the length of fpatch (F1-F2c-F1c/F1-F2-F1c)
F               = 78;
%Specify the length between the patches
D               = 33;

%Enter filename to store the object Arrays (ssobjArray,slobjArray,pdsobjArray,objectArrAccessMat)
mat_fileName_DM = '25SEC_100_F2c_DM.mat';

%-----------------------------------------------------------------------------------------------------------------------------------

%CREATING ARRAY OF OBJECTS
Array       = ones(1,capacity);
ssobjArray  = SSObjectArray(Array);
slobjArray  = SLObjectArray(Array);
pdsobjArray = PDSObjectArray(Array);

%ITERATOR VALUES FOR OBJECT ARRAYS
ssi  = 1;
sli  = 1;
pdsi = 1;

% IMP NOTE : NO OBJECT ARRAY IS CREATED FOR TERMINATED AMPLICONS AS THEY DO NOT PARTICIPATE IN ANY FURTHER REACTIONS 

%THIS IS WHERE PREALLOCATION OF ALL THE ARRAYS AND MATRICES USED IN THE PROGRAM ARE DONE -------------------------------------------

currGen            = 1;
nextGen            = 2;
genNum             = 1;
cg                 = 1;
ng                 = 1;
genArray           = strings(1,colNum);
nextgenArray       = strings(1,colNum);
%formationTimeMat is a matrix which stores the time when the particular amplicon is formed.
formationTimeMat   = zeros(2,colNum); 
%These arrays count the number of pdsa and pdsc amplicons formed respectively
flagpdsaCounterArr = zeros(1,colNum);
flagpdscCounterArr = zeros(1,colNum);
%Stores the index values to access the object arrays for easy validation
objArrAccessMat    = strings(rowNum,colNum);
tIdentifier        = strings(1,colNum);
%Stores the index value to access the amplicons from which the products are formed from
accessArr          = zeros(2,colNum);

numOfRows = ceil(simTime/timeInt) + 1;

%Matrix for storing number of SS,SL,PDS and T formed
%COL 1 : TIME
%COL 2 : NUMBER OF SS FORMED
%COL 3 : NUMBER OF PDS FORMED
%COL 4 : NUMBER OF SL FORMED
%COL 5 : NUMBER OF T FORMED
%COL 6 : NUMBER OF NUCLEOTIDES ADDED
dataMat = zeros(numOfRows,6);

dataMat(1,1) = 0;
dataMat(1,2) = 0;
dataMat(1,3) = 0;
dataMat(1,4) = 0;
dataMat(1,5) = 0;
dataMat(1,6) = 0;

searchTime = timeInt;
for i = 2:numOfRows
    dataMat(i,1) = searchTime;
    searchTime = searchTime + timeInt;
end

%Indices for accessing the object arrays (SS,SL AND PDS)
%In order to create the amplicon we have to get the data from its parent and hence these index values 
%help to access the correct parent amplicon from the object array.

ssAccess  = 1;
slAccess  = 1;
pdsAccess = 1;

%-----------------------------------------------------------------------------------------------------------------------------------



% SECTION 1 : 1st AMPLICON ---------------------------------------------------------------------------------------------------------

%The same variable has been used for iteration in the nucletideAddedMat
formationTimeMat(currGen,cg)   = 0;
genArray(currGen,cg)           = type;

%CHILD
%This is useful to access the SS amplicon in the object array
identifier  = "SS";
arrayIndexC = num2str(ssi);
objArrAccessMat(genNum,cg) = strcat(identifier,arrayIndexC);

numOfCopies = copiesOfDNA;

%This is where the SS amplicon and the amplicons it can form are actually created 
[SSObj,PDSssSS,PDSlsSS,PDSesSS,SLdsSS,naPDSfSS,naSLfSS,naTfSS,numOfPDSfSS] = ssobjFunc(startingSeq,B,F,D);

%GRAND CHILDREN
%Determining what products can be formed and if formed specifying how many and at what time nucleotides are added
formed = 0;
totAmplicons = 2;
slformationTime = formationTimeMat(currGen,cg) + naSLfSS/b;
if slformationTime <= simTime
    formationTimeMat(nextGen,ng)    = slformationTime;
    nextgenArray(1,ng)              = "SLa";
    accessArr(nextGen,ng)           = ssAccess;         
    ng = ng + 1;
    formed = formed + 1;
end

tformationTime  = formationTimeMat(currGen,cg) + naTfSS/b;
if tformationTime <= simTime
    formationTimeMat(nextGen,ng)    = tformationTime;
    nextgenArray(1,ng)              = "T";
    accessArr(nextGen,ng)           = ssAccess;
    ng = ng + 1;
    formed = formed + 1;
end

ssAccess = ssAccess + 1;

%probability of formation of SL is given by
probSL = slformprob;
%probability of formation of other amplicons is given by probOA (probability other amplicons)
probOA = (1 - probSL)/(totAmplicons - 1);

numOfCopiesSL = round(probSL*numOfCopies);
numOfCopiesT  = round(probOA*numOfCopies);

%Specifying end of products formed from SS
formationTimeMat(nextGen,ng)   = 0.01;
nextgenArray(1,ng)             = "|";
ng = ng + 1;

for  i = 2:numOfRows
    if (dataMat(i,1) == round(slformationTime))
        dataMat(i,6) = dataMat(i,6) + naSLfSS*(numOfCopiesSL);
        dataMat(i,4) = dataMat(i,4) + numOfCopiesSL;
    end
    
    if (dataMat(i,1) == round(tformationTime))
        dataMat(i,6) = dataMat(i,6) + naTfSS*(numOfCopiesT); 
        dataMat(i,5) = dataMat(i,5) + numOfCopiesT;
    end
end

%Sending the info to SS object array
%Related to object
ssobjArray(ssi).copies         = numOfCopies;
ssobjArray(ssi).ltSS           = SSObj.tag;
ssobjArray(ssi).ss             = SSObj.ss;
ssobjArray(ssi).PDSssSS        = PDSssSS;
ssobjArray(ssi).PDSlsSS        = PDSlsSS;
ssobjArray(ssi).PDSesSS        = PDSesSS;
ssobjArray(ssi).naPDS          = naPDSfSS;
ssobjArray(ssi).SLdsSS         = SLdsSS;
ssobjArray(ssi).naSL           = naSLfSS;
ssobjArray(ssi).naT            = naTfSS;
ssobjArray(ssi).numOfCopiesPDS = 0;
ssobjArray(ssi).numOfCopiesSL  = numOfCopiesSL;
ssobjArray(ssi).numOfCopiesT   = numOfCopiesT;
ssi = ssi + 1;
    
%genMat stores each generation of the amplicons
genMat = vertcat(genArray,nextgenArray);

totTime  = 0;

firstAmplicon = 1;

%Array index value for lastAmplicons and lastAmpliconsLen arrays
la  = 1;
%These are used to count the copies of SS,SL,PDS and T fomred at the end of the simTime
ss_count  = 0;
sl_count  = 0;
pds_count = 0;
t_count   = numOfCopiesT;

%For resizing purposes
rowBreak = 1;
colBreak = 1;

stopFlag = 0;

% END OF SECTION 1 -----------------------------------------------------------------------------------------------------------------

% MAIN PROBLEM : GETTING THE CORRECT ssAccess, slAccess, pdsAccess. This is
% important so that the next gen amplicons can retreive their data from the
% correct parent.

% didItForm   ->[ '1'    '2'    '0'   '1']
% howManyForm ->[                     '2']
% currentGen  ->[ 'SLa'  'T'    'SSb'  'PDSc'] 
% nextGen     ->[ 'PDSb' 'SS'   'PDSc' 'PDSb']

% In next iteration nextGen -> currentGen and currentGen -> prevGen

% CASE 1 : WHEN THE CURRENT GEN AMPLICON DOES NOT PRODUCE ANY NEXT GEN AMPLICONS
% SOLN   : Use an array called didItForm. This stores only 3 numbers 0,1,2.
% 0 tells that it did not produce any amplicon so when you query the currentGen array and the 
% didItForm array in the nextGen and spot an SS with a value of 0 from the didItForm array then 
% you have to increase ssAccess = ssAccess + 1 so that it can move on to the next SS of the currentGen

% CASE 2 : WHEN THE CURRENT GEN AMPLICON PRODUCES BUT NOT THE FULL SET OF NEXT GEN AMPLICONS
% SOLN   : Take for example an SS which can produce 'PDSa','PDSa','PDSa','SLa','T' (this will be known when SS is created) 
% After this checks are performed to see if all these amplicons are formed within simTime but it happens that only
% 'PDSa','PDSa' are formed. Hence ssAccess = ssAccess + 1 has to be done at the second PDS. In order to know this another 
% array called howManyForm will tell at which child amplicon there has to be an increase in the ssAccess by 1.There wont be cases 
% when PDSc forms 'SS' 'PDSc' (leaves out a 'PDSc') 'SLb'. This wont happen because they are ordered in such a way to reflect 
% increasing time of formation.


% PROBLEM 2 : HOW TO PDSA formed from two different SS ex :SSa & SSb or
% even form same SS like SSa and SSa


% SECTION 2 : MAIN WHILE LOOP ------------------------------------------------------------------------------------------------------

% IMP NOTE : This loop is responsible for creating the amplicons in the current generation 
%            and specifying what all amplicons can be formed in the next generation 

% IMP NOTE : The main loop generates all possible reactions from the amplicons in the previous generation, including the ones which go past the
%            simTime. Hence post processing of the important matrices such as genMat, formationTimeMat, nucleotidesAddedMat has to be done.

while stopFlag == 0
    
tic;

genNum = genNum + 1;

%Assigning the array generated in prev gen to present gen
genArray = nextgenArray;
formationTimeMat(1,:) = formationTimeMat(2,:);
formationTimeMat(2,:) = zeros(1,colNum);
accessArr(1,:)        = accessArr(2,:);
accessArr(2,:)        = zeros(1,colNum);

%some arrays have to be refreshed before they can be used again
nextgenArray       = strings(1,colNum);
pdsaCounterArr     = flagpdsaCounterArr;
pdscCounterArr     = flagpdscCounterArr;
flagpdsaCounterArr = zeros(1,colNum);
flagpdscCounterArr = zeros(1,colNum);

t_Identifier       = tIdentifier;
tIdentifier        = strings(1,colNum);

%Current Gen iterator
cg = 1;
ng = 1;

while true

%----------------------------------------------------------------

    if genArray(1,cg) == "SSa" || genArray(1,cg) == "SSb" || genArray(1,cg) == "SSc"
    
    %PARENT
    %This tells us from which PDS, SS was formed from
    identifierP = "PDS";
    arrayIndexP = num2str(accessArr(currGen,cg));
    fromWhichAmplicon = strcat(identifierP,arrayIndexP);   
        
    %CHILD
    %Creating the SS object
    identifierC = "SS";
    arrayIndexC = num2str(ssi);
    amplicon    = strcat(identifierC,arrayIndexC);
    objArrAccessMat(genNum,cg) = amplicon;
    %This tells us when the SS amplicon was formed
    tagTime = formationTimeMat(currGen,cg); 
    
    %Accesssing the number of copies
    numOfCopies = pdsobjArray(accessArr(currGen,cg)).numOfCopiesSS;
  
    %This is where the SS amplicon and the amplicons it can form are actually created
    [SSObj,PDSssSS,PDSlsSS,PDSesSS,SLdsSS,naPDSfSS,naSLfSS,naTfSS,numOfPDSfSS] = ssobjFunc(pdsobjArray(accessArr(currGen,cg)).ss,B,F,D);
    
    %This shows how many product amplicons can be formed from the parent
    totAmpliconsOA = numOfPDSfSS + 1;
    
    numOfCopiesSL  = round(probSL*numOfCopies);
    numOfCopiesOA  = round((numOfCopies - numOfCopiesSL)/totAmpliconsOA);
    if SSObj.pds == 1
        numOfCopiesPDS = numOfCopiesOA;
    else
        numOfCopiesPDS = 0;
    end
    numOfCopiesT = numOfCopiesOA;
    
    %GRAND CHILDREN
    z = ng;
    formed = 0;
    
    %Specifies at what time PDS will be obtained from SS 
    if SSObj.pds == 1
        for i = 1:numOfPDSfSS
            pdsformationTime(1,i) = formationTimeMat(currGen,cg) + naPDSfSS(1,i)/b;
            totAmplicons = totAmplicons + 1;
            if (pdsformationTime(1,i) <= simTime) && (numOfCopiesPDS > 0) 
                formationTimeMat(nextGen,ng)   = pdsformationTime(1,i);
                nextgenArray(1,ng)             = "PDSa";
                accessArr(nextGen,ng)          = ssAccess;
                ng = ng + 1;
                formed = formed + 1;
            end
        end
    end
    
    %This tells how many PDSA have been formed
    flagpdsaCounterArr(1,z) = formed;
    
    %Specifies at what time SL will be obtained from SS
    slformationTime = formationTimeMat(currGen,cg) + naSLfSS/b;
    totAmplicons = totAmplicons + 1; 
    if (slformationTime <= simTime) && (numOfCopiesSL > 0)
        formationTimeMat(nextGen,ng)   = slformationTime;
        nextgenArray(1,ng)             = "SLa";
        accessArr(nextGen,ng)          = ssAccess;
        ng = ng + 1;
    else
        if numOfCopies > 0
            ss_count = ss_count + numOfCopies;
        end
    end
    
    %Specifies at what time T will be obtained form SS
    tformationTime = formationTimeMat(currGen,cg) + naTfSS/b;
    totAmplicons = totAmplicons + 1;
    if (tformationTime <= simTime) && (numOfCopiesT > 0)
        formationTimeMat(nextGen,ng)   = tformationTime;
        nextgenArray(1,ng)             = "T";
        accessArr(nextGen,ng)          = ssAccess;
        tIdentifier(1,ng)              = strcat("TSS",num2str(ssi));
        ng = ng + 1;
    end
    
    ssAccess = ssAccess + 1;

    cg = cg + 1;
    
    %Specifying end of products formed from SS
    formationTimeMat(nextGen,ng)   = 0.01;
    nextgenArray(1,ng)             = "|";
    ng = ng + 1;
    
    for  i = 2:numOfRows
 
        if SSObj.pds == 1
            for j = 1:length(pdsformationTime)
                if (dataMat(i,1) == round(pdsformationTime(1,j)))
                    dataMat(i,6) = dataMat(i,6) + naPDSfSS(1,j)*(numOfCopiesPDS);
                    dataMat(i,3) = dataMat(i,3) + numOfCopiesPDS;
                end
            end
        end
        
        if (dataMat(i,1) == round(slformationTime))
            dataMat(i,6) = dataMat(i,6) + naSLfSS*(numOfCopiesSL);
            dataMat(i,4) = dataMat(i,4) + numOfCopiesSL; 
        end
        
        if (dataMat(i,1) == round(tformationTime))
            dataMat(i,6) = dataMat(i,6) + naTfSS*(numOfCopiesT);
            dataMat(i,5) = dataMat(i,5) + numOfCopiesT;
        end
        
    end
    
    t_count = t_count + numOfCopiesT;
    
    pdsformationTime = [];
    
    %Sending info about SS to SS object array for storage
    [ssi,ssobjArray] = ssobjArrayFiller(ssobjArray,ssi,SSObj,PDSssSS,PDSlsSS,PDSesSS,SLdsSS,numOfPDSfSS,naTfSS,naSLfSS,naPDSfSS,tagTime,fromWhichAmplicon,numOfCopies,numOfCopiesPDS,numOfCopiesSL,numOfCopiesT);
    
    end
   
%----------------------------------------------------------------

    if genArray(1,cg) == "PDSa" || genArray(1,cg) == "PDSb" || genArray(1,cg) == "PDSc" 
        
        %CREATING THE PDS OBJECT
        if genArray(1,cg) == "PDSa"
            
            x = cg;
            for i = 1:pdsaCounterArr(1,x)
                
                 %PARENT
                 %This tells us from which SS, PDS was formed from
                 identifierP = "SS";
                 arrayIndexP = num2str(accessArr(currGen,cg));
                 fromWhichAmplicon = strcat(identifierP,arrayIndexP);
                
                 %CHILD
                 %Creating the SS object
                 identifierC = "PDS";
                 arrayIndexC = num2str(pdsi);
                 amplicon    = strcat(identifierC,arrayIndexC);
                 objArrAccessMat(genNum,cg) = amplicon;
                 %This tells us when the PDS amplicon was formed
                 tagTime = formationTimeMat(currGen,cg);
                 
                 numOfCopies = ssobjArray(accessArr(currGen,cg)).numOfCopiesPDS;
                
                 %This is where the PDS amplicon and the amplicons it can form are actually created
                 [PDSObj,PDSssPDS,PDSlsPDS,PDSesPDS,SLdsPDS,naPDSfPDS,naSLfPDS,naTfPDS,numOfPDSfPDS] = pdsobjFunc(ssobjArray(accessArr(currGen,cg)).PDSssSS{1,i}, ssobjArray(accessArr(currGen,cg)).PDSlsSS{1,i}, ssobjArray(accessArr(currGen,cg)).PDSesSS,B,F,D);
                 
                 totAmpliconsOA = numOfPDSfPDS + 1;
                 
                  numOfCopiesSS  = 1*numOfCopies;
                  numOfCopiesSL  = round(probSL*numOfCopies);
                  numOfCopiesOA  = round((numOfCopies - numOfCopiesSL)/(totAmpliconsOA));
                  if PDSObj.pds == 1
                      numOfCopiesPDS = numOfCopiesOA;
                  else
                      numOfCopiesPDS = 0;
                  end                     
                  numOfCopiesT  = numOfCopiesOA;
                               
                 %GRAND CHILDREN
                 %Creating the next row of the formationTimeMat
                 [formed,totAmplicons,ssformationTime,pdsformationTime,slformationTime,tformationTime,ng,flagpdscCounterArr,nextgenArray,formationTimeMat,accessArr,pdsAccess,tIdentifierAdd,pds_count] = formationTimeMatGenerator("a",PDSObj,numOfPDSfPDS,simTime,b,naPDSfPDS,naSLfPDS,naTfPDS,cg,ng,flagpdscCounterArr,nextgenArray,formationTimeMat,accessArr,pdsAccess,pds_count,numOfCopies,numOfCopiesSS,numOfCopiesPDS,numOfCopiesSL,numOfCopiesT);
                 
                 if tIdentifierAdd == 1
                     tIdentifier(1,ng - 1) = strcat("TPDS",num2str(pdsi));
                 end
                 
                 cg = cg + 1;
                 
                 %Specifying end of products formed from PDSA
                 formationTimeMat(nextGen,ng)   = 0.01;
                 nextgenArray(1,ng)             = "|";
                 ng = ng + 1;
                 
                 %Sending infor to dataMat
                 for j = 2:numOfRows
                     
                    if (dataMat(j,1) == round(ssformationTime))
                        dataMat(j,2) = dataMat(j,2) + numOfCopiesSS;
                    end
 
                    if pdsformationTime(1,1) ~= 0
                        for k = 1:length(pdsformationTime)
                            if (dataMat(j,1) == round(pdsformationTime(1,k)))
                                dataMat(j,6) = dataMat(j,6) + naPDSfPDS(1,k)*(numOfCopiesPDS);
                                dataMat(j,3) = dataMat(j,3) + numOfCopiesPDS;
                            end
                        end
                    end

                    if (dataMat(j,1) == round(slformationTime))
                        dataMat(j,6) = dataMat(j,6) + naSLfPDS*(numOfCopiesSL);
                        dataMat(j,4) = dataMat(j,4) + numOfCopiesSL;
                    end

                    if (dataMat(j,1) == round(tformationTime))
                        dataMat(j,6) = dataMat(j,6) + naTfPDS*(numOfCopiesT);
                        dataMat(j,5) = dataMat(j,5) + numOfCopiesT;
                    end

                 end
                 
                 tcount = t_count + numOfCopiesT;
                
                 pdsformationTime = [];
                 
                %Sending info to the pds object array for storage
                [pdsi,pdsobjArray] = pdsobjArrayFiller(pdsobjArray,pdsi,PDSObj,PDSssPDS,PDSlsPDS,PDSesPDS,SLdsPDS,numOfPDSfPDS,naTfPDS,naPDSfPDS,naSLfPDS,tagTime,fromWhichAmplicon,numOfCopies,numOfCopiesSS,numOfCopiesPDS,numOfCopiesSL,numOfCopiesT);
                 
            end
            
        elseif genArray(1,cg) == "PDSb"
            
                 %PARENT
                 %This tells us from which SL, PDS was formed from 
                 identifierP = "SL";
                 arrayIndexP = num2str(accessArr(currGen,cg));
                 fromWhichAmplicon = strcat(identifierP,arrayIndexP);
            
                 %CHILD
                 %Creating the PDS object
                 identifierC = "PDS";
                 arrayIndexC = num2str(pdsi);
                 amplicon    = strcat(identifierC,arrayIndexC);
                 objArrAccessMat(genNum,cg) = amplicon;
                 %This tells us when the PDS amplicon was formed
                 tagTime     = formationTimeMat(currGen,cg);
                 
                 numOfCopies = slobjArray(accessArr(currGen,cg)).numOfCopiesPDS;
             
                 %This is where the pds amplicon and the amplicons it can form are actually created
                 [PDSObj,PDSssPDS,PDSlsPDS,PDSesPDS,SLdsPDS,naPDSfPDS,naSLfPDS,naTfPDS,numOfPDSfPDS] = pdsobjFunc(slobjArray(accessArr(currGen,cg)).PDSssSL, slobjArray(accessArr(currGen,cg)).PDSlsSL, slobjArray(accessArr(currGen,cg)).PDSesSL,B,F,D);
                 
                 totAmpliconsOA = numOfPDSfPDS + 1;
                 
                  numOfCopiesSS  = 1*numOfCopies;
                  numOfCopiesSL  = round(probSL*numOfCopies);
                  numOfCopiesOA  = round((numOfCopies - numOfCopiesSL)/(totAmpliconsOA));
                  if PDSObj.pds == 1
                      numOfCopiesPDS = numOfCopiesOA;
                  else
                      numOfCopiesPDS = 0;
                  end                     
                  numOfCopiesT  = numOfCopiesOA;
                 
                 %GRAND CHILDREN
                 %Creating the next row of the formationTimeMat
                 [formed,totAmplicons,ssformationTime,pdsformationTime,slformationTime,tformationTime,ng,flagpdscCounterArr,nextgenArray,formationTimeMat,accessArr,pdsAccess,tIdentifierAdd,pds_count] = formationTimeMatGenerator("b",PDSObj,numOfPDSfPDS,simTime,b,naPDSfPDS,naSLfPDS,naTfPDS,cg,ng,flagpdscCounterArr,nextgenArray,formationTimeMat,accessArr,pdsAccess,pds_count,numOfCopies,numOfCopiesSS,numOfCopiesPDS,numOfCopiesSL,numOfCopiesT);
                 
                 if tIdentifierAdd == 1
                     tIdentifier(1,ng - 1) = strcat("TPDS",num2str(pdsi));
                 end
                 
                 cg = cg + 1;
                 
                 %Specifying end of products formed from PDSA
                 formationTimeMat(nextGen,ng)   = 0.01;
                 nextgenArray(1,ng)             = "|";
                 ng = ng + 1;
                    
                 %Sending info to the dataMat
                 for j = 2:numOfRows
 
                    if (dataMat(j,1) == round(ssformationTime))
                        dataMat(j,2) = dataMat(j,2) + numOfCopiesSS;
                    end
 
                    if pdsformationTime(1,1) ~= 0
                        for k = 1:length(pdsformationTime)
                            if (dataMat(j,1) == round(pdsformationTime(1,k)))
                                dataMat(j,6) = dataMat(j,6) + naPDSfPDS(1,k)*(numOfCopiesPDS);
                                dataMat(j,3) = dataMat(j,3) + numOfCopiesPDS;
                            end
                        end
                    end

                    if (dataMat(j,1) == round(slformationTime))
                        dataMat(j,6) = dataMat(j,6) + naSLfPDS*(numOfCopiesSL);
                        dataMat(j,4) = dataMat(j,4) + numOfCopiesSL;
                    end

                    if (dataMat(j,1) == round(tformationTime))
                        dataMat(j,6) = dataMat(j,6) + naTfPDS*(numOfCopiesT);
                        dataMat(j,5) = dataMat(j,5) + numOfCopiesT;
                    end

                 end
                 
                 pdsformationTime = [];
                 
                 t_count = t_count + numOfCopiesT;
                 
                 %Sending info to the pds object array for storage
                 [pdsi,pdsobjArray] = pdsobjArrayFiller(pdsobjArray,pdsi,PDSObj,PDSssPDS,PDSlsPDS,PDSesPDS,SLdsPDS,numOfPDSfPDS,naTfPDS,naPDSfPDS,naSLfPDS,tagTime,fromWhichAmplicon,numOfCopies,numOfCopiesSS,numOfCopiesPDS,numOfCopiesSL,numOfCopiesT);
                
        else
            
            x = cg;
            for i = 1:pdscCounterArr(1,x)
                
                 %PARENT
                 %This tells us from which PDS, PDS was formed from 
                 identifierP = "PDS";
                 arrayIndexP = num2str(accessArr(currGen,cg));
                 fromWhichAmplicon = strcat(identifierP,arrayIndexP);
                
                 %CHILD
                 %Adding the type of amplicon and its index number so that the amplicon data can be accessed from the pdsobjectArray
                 identifierC = "PDS";
                 arrayIndexC = num2str(pdsi);
                 amplicon    = strcat(identifierC,arrayIndexC);
                 objArrAccessMat(genNum,cg) = amplicon;
                 %This tells us when the PDS amplicon was formed
                 tagTime     = formationTimeMat(currGen,cg);
                 
                 numOfCopies = pdsobjArray(accessArr(currGen,cg)).numOfCopiesPDS;
                 
                 %This is where the pds amplicon and the amplicons it can form are actually created
                 [PDSObj,PDSssPDS,PDSlsPDS,PDSesPDS,SLdsPDS,naPDSfPDS,naSLfPDS,naTfPDS,numOfPDSfPDS] = pdsobjFunc(pdsobjArray(accessArr(currGen,cg)).PDSssPDS{1,i}, pdsobjArray(accessArr(currGen,cg)).PDSlsPDS{1,i}, pdsobjArray(accessArr(currGen,cg)).PDSesPDS,B,F,D);
                 
                 totAmpliconsOA = numOfPDSfPDS + 1;
                 
                  numOfCopiesSS  = 1*numOfCopies;
                  numOfCopiesSL  = round(probSL*numOfCopies);
                  numOfCopiesOA  = round((numOfCopies - numOfCopiesSL)/(totAmpliconsOA));
                  if PDSObj.pds == 1
                      numOfCopiesPDS = numOfCopiesOA;
                  else
                      numOfCopiesPDS = 0;
                  end                     
                  numOfCopiesT  = numOfCopiesOA;

                 %GRAND CHILDREN
                 %Creating the next row of the formationTimeMat
                 [formed,totAmplicons,ssformationTime,pdsformationTime,slformationTime,tformationTime,ng,flagpdscCounterArr,nextgenArray,formationTimeMat,accessArr,pdsAccess,tIdentifierAdd,pds_count] = formationTimeMatGenerator("c",PDSObj,numOfPDSfPDS,simTime,b,naPDSfPDS,naSLfPDS,naTfPDS,cg,ng,flagpdscCounterArr,nextgenArray,formationTimeMat,accessArr,pdsAccess,pds_count,numOfCopies,numOfCopiesSS,numOfCopiesPDS,numOfCopiesSL,numOfCopiesT);
                 
                 if tIdentifierAdd == 1
                     tIdentifier(1,ng - 1) = strcat("TPDS",num2str(pdsi));
                 end
                
                 cg = cg + 1;
                 
                 %Specifying end of products formed from PDSA
                 formationTimeMat(nextGen,ng)   = 0.01;
                 nextgenArray(1,ng)             = "|";
                 ng = ng + 1;
                     
                for j = 2:numOfRows
                    
                    if (dataMat(j,1) == round(ssformationTime))
                        dataMat(j,2) = dataMat(j,2) + numOfCopiesSS;
                    end
 
                    if pdsformationTime(1,1) ~= 0
                        for k = 1:length(pdsformationTime)
                            if (dataMat(j,1) == round(pdsformationTime(1,k)))
                                dataMat(j,6) = dataMat(j,6) + naPDSfPDS(1,k)*(numOfCopiesPDS);
                                dataMat(j,3) = dataMat(j,3) + numOfCopiesPDS;
                            end
                        end
                    end

                    if (dataMat(j,1) == round(slformationTime))
                        dataMat(j,6) = dataMat(j,6) + naSLfPDS*(numOfCopiesSL);
                        dataMat(j,4) = dataMat(j,4) + numOfCopiesSL;
                    end

                    if (dataMat(j,1) == round(tformationTime))
                        dataMat(j,6) = dataMat(j,6) + naTfPDS*(numOfCopiesT);
                        dataMat(j,5) = dataMat(j,5) + numOfCopiesT;
                    end
 
                end
                
                pdsformationTime = [];
                
                t_count = t_count + numOfCopiesT;
                 
                %Sending info to the pds object array for storage
                [pdsi,pdsobjArray] = pdsobjArrayFiller(pdsobjArray,pdsi,PDSObj,PDSssPDS,PDSlsPDS,PDSesPDS,SLdsPDS,numOfPDSfPDS,naTfPDS,naPDSfPDS,naSLfPDS,tagTime,fromWhichAmplicon,numOfCopies,numOfCopiesSS,numOfCopiesPDS,numOfCopiesSL,numOfCopiesT);
                  
            end
            
        end
            
    end

%----------------------------------------------------------------
    
    if genArray(1,cg) == "SLa" || genArray(1,cg) == "SLb"
        
        %CREATING THE SL OBJECT
        if genArray(1,cg) == "SLa"
        
        %PARENT
        identifierP = "SS";
        arrayIndexP = num2str(accessArr(currGen,cg));
        fromWhichAmplicon = strcat(identifierP,arrayIndexP);
        tagTime     = formationTimeMat(currGen,cg);
            
        %CHILD
        identifierC = "SL";
        arrayIndexC = num2str(sli);
        amplicon    = strcat(identifierC,arrayIndexC);
        objArrAccessMat(genNum,cg) = amplicon;
        
        numOfCopies = ssobjArray(accessArr(currGen,cg)).numOfCopiesSL;
       
        [SLObj,PDSssSL,PDSlsSL,PDSesSL,naPDSfSL] = slobjFunc(ssobjArray(accessArr(currGen,cg)).SLdsSS,B,F,D);
       
        else
            
        %PARENT
        identifierP = "PDS";
        arrayIndexP = num2str(accessArr(currGen,cg));
        fromWhichAmplicon = strcat(identifierP,arrayIndexP);
        tagTime     = formationTimeMat(currGen,cg);
        
        %CHILD
        identifierC = "SL";
        arrayIndexC = num2str(sli);
        amplicon    = strcat(identifierC,arrayIndexC);
        objArrAccessMat(genNum,cg) = amplicon;
        
        numOfCopies = pdsobjArray(accessArr(currGen,cg)).numOfCopiesSL;
        
        [SLObj,PDSssSL,PDSlsSL,PDSesSL,naPDSfSL] = slobjFunc(pdsobjArray(accessArr(currGen,cg)).SLdsPDS,B,F,D);
        
        end
        
        if numOfCopies == 0
            numOfCopiesPDS = 0;
        else
            numOfCopiesPDS = numOfCopies;
        end
       
        %Specifies at what time PDS was formed from SL
        pdsformationTime(1,1) = formationTimeMat(currGen,cg) + naPDSfSL/b;
        if (pdsformationTime(1,1) <= simTime) && (numOfCopiesPDS > 0)
            formationTimeMat(nextGen,ng)  = pdsformationTime;
            nextgenArray(1,ng)            = "PDSb";
            accessArr(nextGen,ng)         = slAccess;
            ng = ng + 1;
            formed = formed + 1;
        else
            if numOfCopies > 0
            sl_count = sl_count + numOfCopies;
            end
        end        
        
        slAccess = slAccess + 1;
        
        cg = cg + 1;

        %Specifying end of products formed from SL
        formationTimeMat(nextGen,ng)   = 0.01;
        nextgenArray(1,ng)             = "|";
        ng = ng + 1;

        %Sending info to the dataMat
        for  i = 2:numOfRows
            
            if (dataMat(i,1) == round(pdsformationTime(1,1)))
                dataMat(i,6) = dataMat(i,6) + naPDSfSL*(numOfCopiesPDS);
                dataMat(i,3) = dataMat(i,3) + numOfCopiesPDS;
                break;
            end

        end
        
        pdsformationTime = [];
        
        %Sending info to SL object array for storage
        [sli,slobjArray] = slobjArrayFiller(slobjArray,sli,SLObj,PDSssSL,PDSlsSL,PDSesSL,naPDSfSL,tagTime,fromWhichAmplicon,numOfCopies,numOfCopiesPDS);
        
    end
    
    if genArray(1,cg) == "|"
        objArrAccessMat(genNum,cg) = "|";
        cg = cg + 1;
    end
    
    if genArray(1,cg) == "T"
            if firstAmplicon == 1
                objArrAccessMat(genNum,cg) = "TSS1";
                firstAmplicon = 0; 
            else
                objArrAccessMat(genNum,cg) = t_Identifier(1,cg);
            end
        cg = cg + 1;
    end
    
    if genArray(1,cg) == ""
        break;
    end
    
end

% %Signifies end of array to break out of loop
% howManyFormed(1,cg) = 0.01;

%Required for resizing the matrices
if cg >= colBreak
    colBreak = cg;
end
rowBreak = rowBreak + 1;

%These lines of code stop the main loop
if nextgenArray(1,1) == ""
    stopFlag = 1;
end

% IMP NOTE : THIS IS HOW THE THE MATRIX OF AMPLICONS IS FORMED
if stopFlag == 0
genmat = genMat;
genMat = vertcat(genmat,nextgenArray);
end

disp(".");

toc
tElapsed = toc;
totTime  = totTime + tElapsed;

end

disp(" TOTAL RUNTIME OF THE PROGRAM : ");
disp(totTime);

%Resizing the matrices
reactionNetwork    = genMat(1:rowBreak,1:colBreak);
objArrAccessMat    = objArrAccessMat(1:rowBreak,1:colBreak);

%Adding the nucleotides added before the timeStamp
for v = 2:numOfRows
    dataMat(v,6) = dataMat(v,6) + dataMat(v-1,6);
end

disp(" NUMBER OF SS COPIES AT THE END OF REACTION TIME :")
disp(ss_count)
disp(" NUMBER OF SL COPIES AT THE END OF REACTION TIME :")
disp(sl_count)
disp(" NUMBER OF PDS COPIES AT THE END OF REACTION TIME :") 
disp(pds_count)
disp(" NUMBER OF T COPIES AT THE END OF REACTION TIME :")
disp(t_count)

% END OF SECTION 2 -----------------------------------------------------------------------------------------------------------------



% SAVING IMPORTANT DATA TO BE ACCESSED LATER ---------------------------------------------------------------------------------------

disp("SAVING DATA .........");

save(mat_fileName_DM,'objArrAccessMat','ssobjArray','slobjArray','pdsobjArray','dataMat');
    
fileName = 'reactionNetwork_1.xlsm';

dimensions = [rowBreak,colBreak];
xlswrite(fileName,dimensions,1,'A1:B1');
xlswrite(fileName,reactionNetwork,1,'A2');

% fileName = 'lampData.xls';
% xlswrite(fileName,dataMat,1);


% END OF SECTION 3 -----------------------------------------------------------------------------------------------------------------


% SECTION 4 : FUNCTIONS USED -------------------------------------------------------------------------------------------------------

% --- FUNCTIONS RELATED TO SS ---
function [SSObj,PDSssSS,PDSlsSS,PDSesSS,SLdsSS,naPDSfSS,naSLfSS,naTfSS,numOfPDSfSS] = ssobjFunc(ss,B,F,D)
    
%CREATION OF SS OBJECT
[SSObj,naTfSS,numOfPDSfSS] = SS(ss,B,F,D);

if SSObj.pds == 1
    
    [PDSssSS,PDSlsSS,PDSesSS,naPDSfSS] = SStoPDSGen(SSObj);
    
else
    PDSssSS  = "";
    PDSlsSS  = "";
    PDSesSS  = "";
    naPDSfSS = 0;
end

[SLdsSS,naSLfSS] = SStoSLGen(SSObj);
    
end

function [ssi,ssobjArray] = ssobjArrayFiller(ssobjArray,ssi,SSObj,PDSssSS,PDSlsSS,PDSesSS,SLdsSS,numOfPDSfSS,naTfSS,naSLfSS,naPDSfSS,tagTime,fromWhichAmplicon,numOfCopies,numOfCopiesPDS,numOfCopiesSL,numOfCopiesT)

%SENDING INFO TO OBJECT ARRAYS

% IMP NOTE : PDSssSS AND PDSlsSS ARE CELL ARRAYS WHILE naPDS, textPDS AND timeStampPDS ARE NUMERIC ARRAYS 
        
%RELATED TO SS OBJECT
ssobjArray(ssi).fromWhichAmplicon = fromWhichAmplicon;
ssobjArray(ssi).tagTime           = tagTime;
ssobjArray(ssi).ltSS              = SSObj.tag;
ssobjArray(ssi).ss                = SSObj.ss;
ssobjArray(ssi).copies            = numOfCopies;
%RELATED TO WHAT IT FORMS
ssobjArray(ssi).numOfPDS          = numOfPDSfSS;
ssobjArray(ssi).PDSssSS           = PDSssSS;
ssobjArray(ssi).PDSlsSS           = PDSlsSS;
ssobjArray(ssi).PDSesSS           = PDSesSS;
ssobjArray(ssi).naPDS             = naPDSfSS;
ssobjArray(ssi).SLdsSS            = SLdsSS;
ssobjArray(ssi).naSL              = naSLfSS;
ssobjArray(ssi).naT               = naTfSS;
ssobjArray(ssi).numOfCopiesPDS    = numOfCopiesPDS;
ssobjArray(ssi).numOfCopiesSL     = numOfCopiesSL;
ssobjArray(ssi).numOfCopiesT      = numOfCopiesT;
ssi = ssi + 1;
         
end
% -------------------------------


% --- FUNCTIONS RELATED TO PDS ---
function [PDSObj,PDSssPDS,PDSlsPDS,PDSesPDS,SLdsPDS,naPDSfPDS,naSLfPDS,naTfPDS,numOfPDSfPDS] = pdsobjFunc(ss,ls,es,B,F,D)

%CREATION OF PDS OBJECT
[PDSObj,naTfPDS,numOfPDSfPDS] = PDS(ss,ls,es,B,F,D);

if (PDSObj.pds == 1)
    
    [PDSssPDS,PDSlsPDS,PDSesPDS,naPDSfPDS] = PDStoPDSGen(PDSObj);
    
else
    PDSssPDS   = "";
    PDSlsPDS   = "";
    PDSesPDS   = "";
    naPDSfPDS  = 0;
end

[SLdsPDS,naSLfPDS] = PDStoSLGen(PDSObj);

end

function [formed,totAmplicons,ssformationTime,pdsformationTime,slformationTime,tformationTime,ng,flagpdscCounterArr,nextgenArray,formationTimeMat,accessArr,pdsAccess,tIdentifierAdd,pds_count] = formationTimeMatGenerator(type,PDSObj,numOfPDSfPDS,simTime,b,naPDSfPDS,naSLfPDS,naTfPDS,cg,ng,flagpdscCounterArr,nextgenArray,formationTimeMat,accessArr,pdsAccess,pds_count,numOfCopies,numOfCopiesSS,numOfCopiesPDS,numOfCopiesSL,numOfCopiesT)

formed       = 0;
totAmplicons = 0;

currGen = 1;
nextGen = 2;
%Specifies how many nucleotides and at what time they were added to form SS 
if PDSObj.pds == 1
    ssformationTime = formationTimeMat(currGen,cg) + naPDSfPDS(1,1)/b;
    totAmplicons    = totAmplicons + 1;  
    if (ssformationTime <= simTime) && (numOfCopiesSS > 0)
        formationTimeMat(nextGen,ng) = ssformationTime;
        if type == "a"
            nextgenArray(1,ng) = "SSa";
            accessArr(nextGen,ng)  = pdsAccess;
            ng = ng + 1;
        end
        if type == "b"
            nextgenArray(1,ng) = "SSb";
            accessArr(nextGen,ng)  = pdsAccess;
            ng = ng + 1;
        end
        if type == "c"
            nextgenArray(1,ng) = "SSc";
            accessArr(nextGen,ng)  = pdsAccess;
            ng = ng + 1;
        end
    end
else
    ssformationTime = formationTimeMat(currGen,cg) + naSLfPDS/b;
    totAmplicons    = totAmplicons + 1;
    if (ssformationTime <= simTime) && (numOfCopiesSS > 0)
        formationTimeMat(nextGen,ng)   = ssformationTime;
        if type == "a"
            nextgenArray(1,ng) = "SSa";
            accessArr(nextGen,ng)  = pdsAccess;
            ng = ng + 1;
        end
        if type == "b"
            nextgenArray(1,ng) = "SSb";
            accessArr(nextGen,ng)  = pdsAccess;
            ng = ng + 1;
        end
        if type == "c"
            nextgenArray(1,ng) = "SSc";
            accessArr(nextGen,ng)  = pdsAccess;
            ng = ng + 1;
        end
    end
end

z = ng;
%Specifies how many nucleotides and at what time they were added to form PDS (DOES NOT ALWAYS OCCUR)
if PDSObj.pds == 1
    for i = 1:numOfPDSfPDS
        pdsformationTime(1,i) = formationTimeMat(currGen,cg) + naPDSfPDS(1,i)/b; %#ok<*AGROW>
        totAmplicons = totAmplicons + 1;
        if (pdsformationTime(1,i) <= simTime) && (numOfCopiesPDS > 0)
            formationTimeMat(nextGen,ng)  = pdsformationTime(1,i);
            nextgenArray(1,ng)            = "PDSc";
            accessArr(nextGen,ng)         = pdsAccess;
            ng = ng + 1;
            formed = formed + 1;
        end
    end
else
    pdsformationTime(1,1) = 0;
end

flagpdscCounterArr(1,z) = formed;

%Specifies how many nucleotides and at what time they were added to form SL
slformationTime = formationTimeMat(currGen,cg) + naSLfPDS/b;
totAmplicons    = totAmplicons + 1;
if (slformationTime <= simTime) && (numOfCopiesSL > 0)
    formationTimeMat(nextGen,ng)  = slformationTime;
    nextgenArray(1,ng)            = "SLb";
    accessArr(nextGen,ng)         = pdsAccess;
    ng = ng + 1;
    formed = formed + 1;
else
    if numOfCopies > 0
        pds_count = pds_count + numOfCopies;
    end
end

%Specifies how many nucleotides and at what time they were added to form T
tformationTime = formationTimeMat(currGen,cg) + naTfPDS/b;
totAmplicons   = totAmplicons + 1;
tIdentifierAdd = 0;
if (tformationTime <= simTime) && (numOfCopiesT > 0)
    formationTimeMat(nextGen,ng)  = tformationTime;
    nextgenArray(1,ng)            = "T";
    accessArr(nextGen,ng)         = pdsAccess;
    tIdentifierAdd                = 1;
    ng = ng + 1;
    formed = formed + 1;
end

pdsAccess = pdsAccess + 1;

end

function [pdsi,pdsobjArray] = pdsobjArrayFiller(pdsobjArray,pdsi,PDSObj,PDSssPDS,PDSlsPDS,PDSesPDS,SLdsPDS,numOfPDSfPDS,naTfPDS,naPDSfPDS,naSLfPDS,tagTime,fromWhichAmplicon,numOfCopies,numOfCopiesSS,numOfCopiesPDS,numOfCopiesSL,numOfCopiesT)
                    
%RELATED TO PDS OBJECT
pdsobjArray(pdsi).fromWhichAmplicon = fromWhichAmplicon;
pdsobjArray(pdsi).tagTime           = tagTime;
pdsobjArray(pdsi).ltPDS             = PDSObj.tag;
pdsobjArray(pdsi).nl                = PDSObj.nl;
pdsobjArray(pdsi).ss                = PDSObj.ss;
pdsobjArray(pdsi).ls                = PDSObj.ls;
pdsobjArray(pdsi).es                = PDSObj.es;
pdsobjArray(pdsi).copies            = numOfCopies;
%RELATED TO WHAT IT FORMS
pdsobjArray(pdsi).numOfPDS          = numOfPDSfPDS;
pdsobjArray(pdsi).PDSssPDS          = PDSssPDS;
pdsobjArray(pdsi).PDSlsPDS          = PDSlsPDS;
pdsobjArray(pdsi).PDSesPDS          = PDSesPDS;
pdsobjArray(pdsi).naPDS             = naPDSfPDS;
pdsobjArray(pdsi).SLdsPDS           = SLdsPDS;
pdsobjArray(pdsi).naSL              = naSLfPDS;
pdsobjArray(pdsi).naT               = naTfPDS;
pdsobjArray(pdsi).numOfCopiesSS     = numOfCopiesSS;
pdsobjArray(pdsi).numOfCopiesPDS    = numOfCopiesPDS;
pdsobjArray(pdsi).numOfCopiesSL     = numOfCopiesSL;
pdsobjArray(pdsi).numOfCopiesT      = numOfCopiesT;
pdsi = pdsi + 1;
                    
end
% --------------------------------


% --- FUNCTIONS RELATED TO SL ---
function [SLObj,PDSssSL,PDSlsSL,PDSesSL,naPDSfSL] = slobjFunc(ds,B,F,D)

%CREATION OF SL OBJECT
SLObj = SL(ds,B,F,D);

[PDSssSL,PDSlsSL,PDSesSL,naPDSfSL] = SLtoPDSGen(SLObj);

end

function [sli,slobjArray] = slobjArrayFiller(slobjArray,sli,SLObj,PDSssSL,PDSlsSL,PDSesSL,naPDSfSL,tagTime,fromWhichAmplicon,numOfCopies,numOfCopiesPDS)

%RELATED TO SL OBJECT
slobjArray(sli).fromWhichAmplicon = fromWhichAmplicon;
slobjArray(sli).tagTime           = tagTime;
slobjArray(sli).ltSL              = SLObj.tag;
slobjArray(sli).ds                = SLObj.ds;
slobjArray(sli).copies            = numOfCopies;
%RELATED TO WHAT IT FORMS
slobjArray(sli).PDSssSL           = PDSssSL;
slobjArray(sli).PDSlsSL           = PDSlsSL;
slobjArray(sli).PDSesSL           = PDSesSL;
slobjArray(sli).naPDS             = naPDSfSL;
slobjArray(sli).numOfCopiesPDS    = numOfCopiesPDS;
sli = sli + 1;
        
end
% -------------------------------

     
% END OF SECTION 4 -----------------------------------------------------------------------------------------------------------------
