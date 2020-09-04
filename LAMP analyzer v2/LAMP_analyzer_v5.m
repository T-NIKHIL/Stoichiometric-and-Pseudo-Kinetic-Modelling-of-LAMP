%                             _ _ _ 
% |          /\     |\    /| |     |
% |         /  \    | \  / | |     |
% |        /_ _ \   |  \/  | |_ _ _|
% |       /      \  |      | |
% |_ _ _ /        \ |      | |

% THIS IS A PROGRAM THAT ANALYZES THE EXPONENTIAL PHASE OF LAMP 

% NAME OF PROGRAM : LAMP_Analyzer_v2

% This is version 2.0 of LAMP_Analyzer and uses parallelization techniques to obtain data for simulation times graeter than 1 minute

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

% SIM TIME - REAL TIME(MIN)
%     1    -    1
%     2    -    7
%     5    -  

%USER INPUT 
inputSeq    = ["B2","F2c"];
compSeq     = ["F2","B2c"];
%Specify in the order as [B;F;D]
patchLen    = [95;78;33];
%Specify the initial number of 'B' patches and 'F' patches
bCount      = 1;
fCount      = 1;
intPatch    = [bCount,fCount,(bCount + fCount - 1)];
%Specify the rate of nucleotide incorporation (nucleotides/second) by Polymerase
rate        = 120;
simTime     = 60;
% DO NOT CHANGE timeInt !
timeInt     = 1;
numOfRows   = simTime;
rowsDataMat = numOfRows; 
%Specify the probability of formation of SL and copies of DNA
probSL      = 0.67;
%Specify the initial number of copies of DNA
copiesOfDNA = 100;
numOfRowsBuffer = 100;
replenishRowAmt = 100;
%Specify if you want to double the number of functions running parallely
increaseNumFun  = 0;
%Specify the file name to save the .mat file
specifyFileName = '1MIN100COPIES_F2c';
% ALLOCATING SPACE TO ARRAYS

% LIST OF IMPORTANT ARRAYS :
% 1) formationTimeArr      : stores the formation times of the products belonging to the dumbell loop
% 2) nucleotidesAddedArr   : stores the nucleotides added to form the products belonging to the dumbell loop
% 3) lastAmpliconsLen      : stores the length of terminated amplicons and the amplicons formed at the end of simulation time
% 4) dataMat               : data repository

%Matrix fo storing number of SS,SL,PDS and T formed
%COL 1 : TIME
%COL 2 : NUMBER OF SS FORMED
%COL 3 : NUMBER OF PDS FORMED
%COL 4 : NUMBER OF SL FORMED
%COL 5 : NUMBER OF T FORMED
%COL 6 : NUMBER OF NUCLEOTIDES ADDED
dataMat             = zeros(numOfRows,6);
dataMatBlank        = zeros(numOfRows,6);
formationTimeArr    = zeros(1,6);
la = 1;

dataMat(1,1)      = 0;
dataMatBlank(1,1) = 0; 
dataMat(1,2)      = 0;
dataMatBlank(1,2) = 0;
dataMat(1,3)      = 0;
dataMatBlank(1,3) = 0;
dataMat(1,4)      = 0;
dataMatBlank(1,4) = 0;
dataMat(1,5)      = 0;
dataMatBlank(1,5) = 0;
dataMat(1,6)      = 0;
dataMatBlank(1,6) = 0;

% NOTE : lastAmpliconsLen stores reactants which have
% 1) formed at the reaction time
% 2) formed before the reaction time but could not form the products within
% the simulation time

% This allocates the time interval values into the time column of the dataMat
searchTime = timeInt;
for i = 2:numOfRows
    dataMat(i,1)      = searchTime;
    dataMatBlank(i,1) = searchTime; 
    searchTime        = searchTime + timeInt;
end

% Variables used to count the number of SS,SL,PDS and T amplicons produced at the end of simTime
ss_copies_count  = 0;
sl_copies_count  = 0;
pds_copies_count = 0;
t_copies_count   = 0;

% NOTES : 
% 1) LENGTH OF AMPLICONS
% 1.1 Length of SL amplicon  = length of es sequence of PDS
% 1.2 Length of PDS amplicon = length of ss sequence of PDS (the double stranded region of PDS)
% 1.3 Length of SS amplicon  = length of ss sequence of SS
% 1.4 Length of T amplicon   = length of ss sequence (if formed from SS) OR length of es sequence (if formed from PDS)
% 2) LENGTHS OF AMPLICONS ARE CALCULATED BY USING THE NUMBER OF 'B' PATCHES AND 'F' PATCHES PRESENT IN THE AMPLICON
% 3) SEQUENCES ARE ONLY USED FOR FINDING THE PRIMER BINDING SITES. ALL SEQUENCES HAVE THEIR 3' AT THE RIGHT END 

% DUMBELL LOOP 
% SS(ss)->SL->PDS->SS(cs)->SL->PDS->SS(ss)

% PRODUCTS GENERATED FROM THIS LOOP WHICH ARE NOT PARTICIPATING IN THE LOOP ARE 
% T from SS, T from PDS, SL from PDS

% Two important arrays are slTimeArr1,slTimeArr2,slCopiesArr1,slCopiesArr2
ta1 = 1;
ca1 = 1;
ta2 = 1;
ca2 = 1;

lengthSS  = intPatch*patchLen;
lengthSL  = intPatch*patchLen;
% nucleotides added to form PDS
naPDS     = lengthSS;

% 1. SS
ss  = inputSeq;
css = compSeq;
formationTimeArr(1,1)  = 0;
extensionTimeArr(1,1)  = 0;

% 2. SS->SL & SS->T
% Here we calculate the time taken to form SL and the number of nucleotides used to form T
if inputSeq(1,end) == "F2c"
naSL                      = [bCount,fCount - 1,(bCount + fCount - 1)]*patchLen;
else
naSL                      = [bCount - 1,fCount,(bCount + fCount - 1)]*patchLen;
end
naTfSS                    = [bCount,fCount,(bCount + fCount - 1)]*patchLen;
extensionTimeArr(1,2)     = (naSL/rate);
formationTimeArr(1,2)     = extensionTimeArr(1,2) + formationTimeArr(1,1);

% 3. SL->PDS
% Here we calculate the time taken to form PDS
extensionTimeArr(1,3)     = (naPDS/rate);
formationTimeArr(1,3)     = extensionTimeArr(1,3) + formationTimeArr(1,2);

% 4. PDS->SS,PDS->SL & PDS->T
% Here we calculate the time taken to form the complementary strand SS and also the number of nucleotides used to form T from PDS
ssPDS                     = css;
cssPDS                    = ss;
esPDS                     = [cssPDS,ssPDS(2:end)];
% Creating the structures of SL 
ds1                       = [esPDS;fliplr(ssPDS),""];
naSLfPDS1                 = [bCount,fCount,(bCount + fCount)]*patchLen;
if inputSeq(1,end) == "F2c"
naTfPDS1                  = [2*bCount,fCount,(2*bCount + fCount - 1)]*patchLen;
else
naTfPDS1                  = [bCount,2*fCount,(bCount + 2*fCount - 1)]*patchLen;
end
lengthPDS1                = naTfPDS1;
extensionTimeArr(1,4)     = (naSLfPDS1/rate);
formationTimeArr(1,4)     = extensionTimeArr(1,4) + formationTimeArr(1,3);


% 5. SS->SL & SS->T
%Here we calculate the time taken to form SL
%Number of nucleotides added to form T from this complementary SS is the same as the initial SS
if inputSeq(1,end) == "F2c"
naSL                      = [bCount - 1,fCount,(bCount + fCount - 1)]*patchLen;
else
naSL                      = [bCount,fCount - 1,(bCount + fCount - 1)]*patchLen;
end
extensionTimeArr(1,5)     = (naSL/rate);
formationTimeArr(1,5)     = extensionTimeArr(1,5) + formationTimeArr(1,4);

% 6. SL->PDS
%Here we calculate the time taken to form PDS
extensionTimeArr(1,6)     = (naPDS/rate);
formationTimeArr(1,6)     = extensionTimeArr(1,6) + formationTimeArr(1,5);

% 7. PDS->SS,PDS->SL & PDS->T 
%Here we calculate the time taken to form SS back again and also the number of nucleotides used to form T from PDS
ssPDS                     = ss;
cssPDS                    = css;
esPDS                     = [cssPDS,ssPDS(2:end)];
ds2                       = [esPDS;fliplr(ssPDS),""];
naSLfPDS2                 = [bCount,fCount,(bCount + fCount)]*patchLen;
if inputSeq(1,end) == "F2c"
naTfPDS2                  = [bCount,2*fCount,(bCount + 2*fCount - 1)]*patchLen;
else
naTfPDS2                  = [2*bCount,fCount,(2*bCount + fCount - 1)]*patchLen;
end
lengthPDS2                = naTfPDS2;
extensionTimeArr(1,1)     = (naSLfPDS2/rate);


% NOTE :
% 1) dataMat stores the numbers of copies of each type of product formed at a particular time
% 2) lastAmpliconsLen stores only the copies of those products formed just before simulation time runs out and also the reactants that are still
%    undergoing transformation to form products 

% SENDING INFO INTO dataMat

% NOTE:
% 1. WHEN PROBABILITY OF FORMATION OF SL IS HIGH (GREATER THAN OR EQUAL TO 50%)
% 1.1 Then a stage will reach when the number of copies of SS,SL and PDS do not change and remain constant. This occurs because of the rounding effect 
% and when you have low number of copies of SS and these produce SL entirely and all these copies of SL will form PDS. All these copies of PDS will inturn
% produce SL and so on. In this case the loop terminates because one of the amplicons participating in the loop has a formation time exceeding that of the simulation time.
% 1 SS -> (probSL*1) = 1 SL -> 1 PDS -> 1 SS -> ....
%
% 2. WHEN PROBABILITY OF FORMATION OF SL IS LOW (LESSER THAN 50%)
% 2.1 In this case the loop terminates because one of the amplicons participating in the loop has reached 0 number of copies. The model will fail when probSL is given 
% values of 50% or lower. This because checks for number of copies for SS,SL or PDS is not done. The reason why its not done is because in reality probSL will always be above 80%

first = 1;

 while true
    
if first == 1
    numOfCopiesSS = copiesOfDNA;
    first = 0;
end

%SENDING SS INFO (PDS->SS)
ssformationTime = round(formationTimeArr(1,1));
index = ssformationTime/timeInt + 1;
if index <= simTime
   dataMat(index,2) = dataMat(index,2) + numOfCopiesSS;
end
if formationTimeArr(1,1) > simTime
    pds_copies_count = pds_copies_count + numOfCopiesPDS;
    lastAmpliconsLen(1,la) = lengthPDS2; %#ok<*SAGROW>
    lastAmpliconsLen(2,la) = numOfCopiesPDS;
    la = la + 1;
    break;
end

% -------------------------------------------
numOfCopiesSL = round(probSL*numOfCopiesSS);
numOfCopiesT  = numOfCopiesSS - numOfCopiesSL;
% -------------------------------------------
    
%SENDING SL INFO (SS->SL)
slformationTime = round(formationTimeArr(1,2));
index = slformationTime/timeInt + 1;
if index <= simTime
    dataMat(index,4) = dataMat(index,4) + numOfCopiesSL;
    dataMat(index,6) = dataMat(index,6) + extensionTimeArr(1,2)*120*numOfCopiesSL;    
end
if formationTimeArr(1,2) > simTime
    ss_copies_count = ss_copies_count + numOfCopiesSS;
    lastAmpliconsLen(1,la) = lengthSS;
    lastAmpliconsLen(2,la) = numOfCopiesSS;
    la = la + 1;
    break; 
end

%NOTE : If SL cannot be formed within the simulation time then T cannot be formed within the simulation time
%SENDING T INFO (SS->T) (TERMINATION PRODUCT)
tformationTime = round(formationTimeArr(1,1) + naTfSS/rate);
index = tformationTime/timeInt + 1;
if index <= simTime
    dataMat(index,5) = dataMat(index,5) + numOfCopiesT;
    dataMat(index,6) = dataMat(index,6) + naTfSS*numOfCopiesT;
end
if ((formationTimeArr(1,1) + naTfSS/rate) <= simTime) && (numOfCopiesT > 0) 
    t_copies_count = t_copies_count + numOfCopiesT;
    lastAmpliconsLen(1,la) = naTfSS;
    lastAmpliconsLen(2,la) = numOfCopiesT;
    la = la + 1;
end

% -------------------------------------------
numOfCopiesPDS = numOfCopiesSL;
% -------------------------------------------

%SENDING PDS INFO (SL->PDS)
pdsformationTime = round(formationTimeArr(1,3));
index = pdsformationTime/timeInt + 1;
if index <= simTime
    dataMat(index,3) = dataMat(index,3) + numOfCopiesPDS;
    dataMat(index,6) = dataMat(index,6) + extensionTimeArr(1,3)*120*numOfCopiesPDS;    
end
if formationTimeArr(1,3) > simTime
    sl_copies_count = sl_copies_count + numOfCopiesSL;
    lastAmpliconsLen(1,la) = lengthSL;
    lastAmpliconsLen(2,la) = numOfCopiesSL; 
    la = la + 1;
    break;  
end

% -------------------------------------------
numOfCopiesSS  = numOfCopiesPDS;
numOfCopiesSL  = round(probSL*numOfCopiesPDS);
numOfCopiesT   = numOfCopiesPDS - numOfCopiesSL;
% -------------------------------------------

% SENDING SL INFO INTO DATA CONTAINER
slTimeArr1(ta1,1)   = (formationTimeArr(1,3) + (naSLfPDS1/rate));
slCopiesArr1(ca1,1) = numOfCopiesSL;
ta1 = ta1 + 1;
ca1 = ca1 + 1;

% SENDING SL INFO INTO dataMat
slformationTime = round(formationTimeArr(1,3) + extensionTimeArr(1,4));
index = slformationTime/timeInt + 1;
if (index <= simTime)
    dataMat(index,4) = dataMat(index,4) + numOfCopiesSL;
    dataMat(index,6) = dataMat(index,6) + naSLfPDS1*numOfCopiesSL;
end

% SENDING SS INFO (PDS->SS)
ssformationTime = round(formationTimeArr(1,4));
index = ssformationTime/timeInt + 1;
if index <= simTime
    dataMat(index,2) = dataMat(index,2) + numOfCopiesSS;
end
if (formationTimeArr(1,4) > simTime)
    pds_copies_count = pds_copies_count + numOfCopiesPDS;
    lastAmpliconsLen(1,la) = lengthPDS1;
    lastAmpliconsLen(2,la) = numOfCopiesPDS;
    la = la + 1;
    break;
end

% NOTE : If SS cannot be formed within the simulation time then T cannot be formed within the simulation time
% SENDING T INFO (PDS->T) (TERMINATION PRODUCT)
tformationTime = round(formationTimeArr(1,3) + naTfPDS1/rate);
index = tformationTime/timeInt + 1;
if index <= simTime
    dataMat(index,5) = dataMat(index,5) + numOfCopiesT;
    dataMat(index,6) = dataMat(index,6) + naTfPDS1*numOfCopiesT;
end
if (((formationTimeArr(1,3) + naTfPDS1/rate) <= simTime) && (numOfCopiesT > 0))
    t_copies_count = t_copies_count + numOfCopiesT;
    lastAmpliconsLen(1,la) = naTfPDS1;
    lastAmpliconsLen(2,la) = numOfCopiesT;
    la = la + 1;
end

% ------------------------------------------
numOfCopiesSL  = round(probSL*numOfCopiesSS);
numOfCopiesT   = numOfCopiesSS - numOfCopiesSL;
% ------------------------------------------

%SENDING SL INFO (SS->SL)
slformationTime = round(formationTimeArr(1,5));
index = slformationTime/timeInt + 1;
if index <= simTime
    dataMat(index,4) = dataMat(index,4) + numOfCopiesSL;
    dataMat(index,6) = dataMat(index,6) + extensionTimeArr(1,5)*120*numOfCopiesSL;
end
if formationTimeArr(1,5) > simTime
    ss_copies_count = ss_copies_count + numOfCopiesSS;
    lastAmpliconsLen(1,la) = lengthSS;
    lastAmpliconsLen(2,la) = numOfCopiesSS;
    la = la + 1;
    break;
end

%NOTE : If SL cannot be formed within the simulation time then T cannot be formed within the simulation time
%SENDING T INFO (SS->T) (TERMINATION PRODUCT)
tformationTime = round(formationTimeArr(1,4) + naTfSS/rate);
index = tformationTime/timeInt + 1;
if index <= simTime
    dataMat(index,5) = dataMat(index,5) + numOfCopiesT;
    dataMat(index,6) = dataMat(index,6) + naTfSS*numOfCopiesT;
end
if ((formationTimeArr(1,4) + naTfSS/rate) <= simTime) && (numOfCopiesT > 0)
    t_copies_count = t_copies_count + numOfCopiesT;
    lastAmpliconsLen(1,la) = naTfSS;
    lastAmpliconsLen(2,la) = numOfCopiesT;
    la = la + 1;
end

% ------------------------------------------
numOfCopiesPDS = numOfCopiesSL;
% ------------------------------------------

%SENDING PDS INFO (SL->PDS)

pdsformationTime = round(formationTimeArr(1,6));
index = pdsformationTime/timeInt + 1;
if index <= simTime
    dataMat(index,3) = dataMat(index,3) + numOfCopiesPDS;
    dataMat(index,6) = dataMat(index,6) + extensionTimeArr(1,6)*120*numOfCopiesPDS;      
end
if formationTimeArr(1,6) > simTime
    sl_copies_count = sl_copies_count + numOfCopiesSL;
    lastAmpliconsLen(1,la) = lengthSL;
    lastAmpliconsLen(2,la) = numOfCopiesSL;
    la = la + 1;
    break;
end

% ------------------------------------------
numOfCopiesSS  = numOfCopiesPDS; 
numOfCopiesSL  = round(probSL*numOfCopiesPDS);
numOfCopiesT   = numOfCopiesPDS - numOfCopiesSL;
% ------------------------------------------

% Sending SL info into data container
slTimeArr2(ta2,1)   = (formationTimeArr(1,6) + (naSLfPDS2/rate));
slCopiesArr2(ca2,1) = numOfCopiesSL;
ta2 = ta2 + 1;
ca2 = ca2 + 1;

% SENDING SL INFO (PDS->SL) 
slformationTime = round(formationTimeArr(1,6) + (naSLfPDS2/rate));
index = slformationTime/timeInt + 1;
if (index <= simTime)
    dataMat(index,4) = dataMat(index,4) + numOfCopiesSL;
    dataMat(index,6) = dataMat(index,6) + naSLfPDS2*numOfCopiesSL;                                        
end

% SENDING T INFO (PDS->T) (TERMINATION PRODUCT)
tformationTime = round(formationTimeArr(1,6) + (naTfPDS2/rate));
index = tformationTime/timeInt + 1;
if index <= simTime
    dataMat(index,5) = dataMat(index,5) + numOfCopiesT;
    dataMat(index,6) = dataMat(index,6) + naTfPDS2*numOfCopiesT;
end
if ((formationTimeArr(1,6) + (naTfPDS2/rate)) <= simTime) && (numOfCopiesT > 0)
    t_copies_count = t_copies_count + numOfCopiesT;
    lastAmpliconsLen(1,la) = naTfPDS2;
    lastAmpliconsLen(2,la) = numOfCopiesT;
    la = la + 1;
end

formationTimeArr(1,1) = formationTimeArr(1,6) + extensionTimeArr(1,1);
formationTimeArr(1,2) = formationTimeArr(1,1) + extensionTimeArr(1,2);
formationTimeArr(1,3) = formationTimeArr(1,2) + extensionTimeArr(1,3);
formationTimeArr(1,4) = formationTimeArr(1,3) + extensionTimeArr(1,4);
formationTimeArr(1,5) = formationTimeArr(1,4) + extensionTimeArr(1,5);
formationTimeArr(1,6) = formationTimeArr(1,5) + extensionTimeArr(1,6);

 end
 
 

% ----------------------------------- DECLARING VARIABLES --------------------------------

% Creating the package buffers (these will store the SS,PDS and SL packages which will then be sent into the respective functions
ssPackageBuffer   = cell(numOfRowsBuffer,2);
slPackageBuffer2  = cell(numOfRowsBuffer,2);
pdsPackageBuffer  = cell(numOfRowsBuffer,2);
slPackageBuffer1  = cell(numOfRowsBuffer,2);

% These are the indices used in the packageBuffers
% SL1 and SL2 are used to differentiate between the two SL functions
% P1  and P2 represent the the different pathways coming out from the dumbell loop
P1bufferSL1 = 1;
P2bufferSL1 = 1;
P1bufferSL2 = 1;
P2bufferSL2 = 1;
P1bufferSS  = 1;
P2bufferSS  = 1;
P1bufferPDS = 1;
P2bufferPDS = 1;

% ----------------------------------------------------------------------------------------

% Sending the initial data first to sltopdsHwy
slPackage1{1} = ds1;
slPackage1{2} = slTimeArr1;
slPackage1{3} = slCopiesArr1;
slPackage2{1} = ds2;
slPackage2{2} = slTimeArr2;
slPackage2{3} = slCopiesArr2;

dataFun = slHwy(slPackage1,slPackage2,dataMatBlank,simTime,timeInt,rate,probSL,patchLen);

% STEP 1 : OBTAINING THE LENGTH OF dataFun
if length(dataFun) == 7
    numAmpliconInput = 1;
    lastData         = dataFun{6};
    addData          = dataFun{7};  
else
    numAmpliconInput = 2;
    lastData         = dataFun{13};
    addData          = dataFun{14};
end

% STEP 2 : ADDING THE addData OBTAINED FROM dataFun INTO COMMON dataMat
for x = 2:6
    dataMat(1:end,x) = dataMat(1:end,x) + addData(1:end,x);
end

% STEP 3 : ADDING THE lastData OBTAINED FROM dataFun INTO COMMON dataMat
if isempty(lastData) == 0
    
    dim = size(lastData);
    lenOfArr = dim(1,2);
    
    for i = 1:lenOfArr
        if lastData(1,i) == 1
            ss_copies_count  = ss_copies_count  + lastData(3,i);
        elseif lastData(1,i) == 2
            pds_copies_count = pds_copies_count + lastData(3,i); 
        elseif lastData(1,i) == 3
            sl_copies_count  = sl_copies_count  + lastData(3,i);
        else
            t_copies_count   = t_copies_count   + lastData(3,i);
        end
    end
    
    lastAmpliconsLen = [lastAmpliconsLen,lastData(2:3,:)];
    lengths          = lastAmpliconsLen';
    lengths          = sortrows(lengths,1);
    [r,c]            = size(lengths);

    uniqueLens       = unique(lengths(:,1));
    countOfLen       = hist(lengths(:,1),uniqueLens);

    len = 1;
    m   = 1;

    while true
        lengthData(m,1) = lengths(len,1);
        lengthData(m,2) = sum(lengths(len:(len + countOfLen(1,m) - 1),2));
        len = len + countOfLen(1,m);
        m  = m + 1; 
        if len > r
            break;
        end
    end
    lastAmpliconsLen = lengthData';
end

% STEP 4 : EXTRACTING THE SS AND PDS DATA
for x = 1:numAmpliconInput

% First unpacking the SS from PDS data
if x == 1
    strSSfPDS   = dataFun{1};
    dataSSfPDS  = dataFun{2};
else
    strSSfPDS   = dataFun{8};
    dataSSfPDS  = dataFun{9};
end

if isempty(strSSfPDS) == 0
    dim         = size(dataSSfPDS);
    numOfCol    = dim(1,2);
    numOfIter   = numOfCol/2;
    shiftCol    = 0;

    % Creating the SS buffer array
    for y = 1:numOfIter
        if x == 1
            ssPackage1{1}        = strSSfPDS{1,y};
            ssPackage1{2}        = strSSfPDS{2,y};
            ssTimeArr            = dataSSfPDS(1:end,1 + shiftCol);
            ssPackage1{3}        = ssTimeArr;
            ssPackage1{4}        = dataSSfPDS(1:end,2 + shiftCol);
            if ssTimeArr(1,1) <= simTime
                ssPackageBuffer{P1bufferSS,1} = ssPackage1;
                P1bufferSS = P1bufferSS + 1;
            end
        else
            ssPackage2{1}        = strSSfPDS{1,y};
            ssPackage2{2}        = strSSfPDS{2,y};
            ssTimeArr            = dataSSfPDS(1:end,1 + shiftCol);
            ssPackage2{3}        = ssTimeArr;
            ssPackage2{4}        = dataSSfPDS(1:end,2 + shiftCol);
            if ssTimeArr(1,1) <= simTime
                ssPackageBuffer{P2bufferSS,2} = ssPackage2;
                P2bufferSS = P2bufferSS + 1;
            end
        end
        shiftCol = shiftCol + 2;
    end
end

% Lastly unpacking the PDS from PDS data
if x == 1
    strPDSfPDS  = dataFun{3};
    dataPDSfPDS = dataFun{4};
    childPDSnum = dataFun{5};
else
    strPDSfPDS  = dataFun{10};
    dataPDSfPDS = dataFun{11};
    childPDSnum = dataFun{12};
end

lengthOfArr = length(childPDSnum);
shiftCol1   = 0;
shiftCol2   = 0;

if isempty(strPDSfPDS) == 0
    % Creating the PDS buffer array
    for y = 1:lengthOfArr
        numOfChildPDS   = childPDSnum(1,y);
        pdsCopiesArr    = dataPDSfPDS(1:end,(numOfChildPDS + 1 + shiftCol2));
        for z = 1:numOfChildPDS
            if x == 1
                pdsPackage1{1}          = strPDSfPDS{1,z + shiftCol1};
                pdsPackage1{2}          = strPDSfPDS{2,z + shiftCol1};
                pdsTimeArr              = dataPDSfPDS(1:end,z + shiftCol2);
                pdsPackage1{3}          = pdsTimeArr;
                pdsPackage1{4}          = pdsCopiesArr;
                if pdsTimeArr(1,1) <= simTime
                    pdsPackageBuffer{P1bufferPDS,1} = pdsPackage1;
                    P1bufferPDS = P1bufferPDS + 1;
                end
            else
                pdsPackage2{1}          = strPDSfPDS{1,z + shiftCol1};
                pdsPackage2{2}          = strPDSfPDS{2,z + shiftCol1};
                pdsTimeArr              = dataPDSfPDS(1:end,z + shiftCol2);
                pdsPackage2{3}          = pdsTimeArr;
                pdsPackage2{4}          = pdsCopiesArr;
                if pdsTimeArr(1,1) <= simTime
                    pdsPackageBuffer{P2bufferPDS,2} = pdsPackage2;
                    P2bufferPDS = P2bufferPDS + 1;
                end
            end
        end
        shiftCol1 = shiftCol1 + numOfChildPDS;
        shiftCol2 = shiftCol2 + numOfChildPDS + 1;
    end
end

end

% ----------------------------------- PARENT WHILE LOOP ----------------------------------

while true
% --------------------------------------- SECTION 1 --------------------------------------

% This gets activated when space in the ssPackageBuffer is running low
if increaseNumFun == 1
    numPackages = 40;
else
    numPackages = 8;
end
dim = size(ssPackageBuffer);
if dim(1,1) <= numPackages
    transferCell    = ssPackageBuffer(1:end,:); 
    ssPackageBuffer = cell(replenishRowAmt,2);
    ssPackageBuffer(1:dim(1,1),:) = transferCell;
end

% Finding how many packages are left in the SS buffer
findEmptySlot   = find(cellfun(@isempty,ssPackageBuffer),1);
lastRowSS       = findEmptySlot - 1;
if isempty(lastRowSS) == 1
    lastRowSS  = dim(1,1); 
end

if lastRowSS ~= 0

if increaseNumFun == 1
    
    % Running 40 functions simultaneuosly
    funcs     = {@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop};
    arguments = {ssPackageBuffer{1,1},ssPackageBuffer{1,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{2,1},ssPackageBuffer{2,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{3,1},ssPackageBuffer{3,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{4,1},ssPackageBuffer{4,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{5,1},ssPackageBuffer{5,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{6,1},ssPackageBuffer{6,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{7,1},ssPackageBuffer{7,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{8,1},ssPackageBuffer{8,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{9,1},ssPackageBuffer{9,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{10,1},ssPackageBuffer{10,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{11,1},ssPackageBuffer{11,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{12,1},ssPackageBuffer{12,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{13,1},ssPackageBuffer{13,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{14,1},ssPackageBuffer{14,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{15,1},ssPackageBuffer{15,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{16,1},ssPackageBuffer{16,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{17,1},ssPackageBuffer{17,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{18,1},ssPackageBuffer{18,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{19,1},ssPackageBuffer{19,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{20,1},ssPackageBuffer{20,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{21,1},ssPackageBuffer{21,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{22,1},ssPackageBuffer{22,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{23,1},ssPackageBuffer{23,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{24,1},ssPackageBuffer{24,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{25,1},ssPackageBuffer{25,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{26,1},ssPackageBuffer{26,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{27,1},ssPackageBuffer{27,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{28,1},ssPackageBuffer{28,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{29,1},ssPackageBuffer{29,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{30,1},ssPackageBuffer{30,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{31,1},ssPackageBuffer{31,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{32,1},ssPackageBuffer{32,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{33,1},ssPackageBuffer{33,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{34,1},ssPackageBuffer{34,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{35,1},ssPackageBuffer{35,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{36,1},ssPackageBuffer{36,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{37,1},ssPackageBuffer{37,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{38,1},ssPackageBuffer{38,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{39,1},ssPackageBuffer{39,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{40,1},ssPackageBuffer{40,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen};
    solutions = cell(40,1); 
    
    % Deleting the packages that have been sent into the parfor to clear out space in the RAM
    ssPackageBuffer(1:40,:) = [];
    P1bufferSS = P1bufferSS - 40;
    P2bufferSS = P2bufferSS - 40;
    if P1bufferSS <= 0
        P1bufferSS = 1;
    end
    if P2bufferSS <= 0
        P2bufferSS = 1;
    end
    
    % Here the 40 functions run parallely
    parfor ii = 1:40
        solutions{ii} = funcs{ii}(arguments{ii,:});
    end
    
    numOfIterX = 40;
    
else
    
    % Running 8 functions simultaneuosly
    funcs     = {@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop,@ssLoop};
    arguments = {ssPackageBuffer{1,1},ssPackageBuffer{1,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{2,1},ssPackageBuffer{2,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{3,1},ssPackageBuffer{3,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{4,1},ssPackageBuffer{4,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{5,1},ssPackageBuffer{5,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{6,1},ssPackageBuffer{6,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{7,1},ssPackageBuffer{7,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 ssPackageBuffer{8,1},ssPackageBuffer{8,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen};
    solutions = cell(8,1); 
    
    % Deleting the packages that have been sent into the parfor
    ssPackageBuffer(1:8,:)  = [];
    P1bufferSS = P1bufferSS - 8;
    P2bufferSS = P2bufferSS - 8;
    if P1bufferSS <= 0
        P1bufferSS = 1;
    end
    if P2bufferSS <= 0
        P2bufferSS = 1;
    end
    
    % Here the 8 functions run parallely
    parfor ii = 1:8
        solutions{ii} = funcs{ii}(arguments{ii,:});
    end
    
    numOfIterX = 8;
    
end

% STEP 1 : EXTRACTING THE DATA FROM SOLUTIONS
for x = 1:numOfIterX
    
    dataFun = solutions{x};
   
    % STEP 2 : DETERMINING THE LENGTH OF dataFun
    if length(dataFun) == 15
        numAmpliconInput = 1;
        lastData         = dataFun{14};
        addData          = dataFun{15};
    elseif length(dataFun) == 30
        numAmpliconInput = 2;
        lastData         = dataFun{29};
        addData          = dataFun{30};
    else
        continue;
    end
    
    % STEP 3 : ADDING THE addData INTO dataMat
    for y = 2:6
        dataMat(1:end,y) = dataMat(1:end,y) + addData(1:end,y);
    end
    
    % STEP 4 : ADDING THE lastData INTO lastAmpliconsLen
    if isempty(lastData) == 0
        
        dim = size(lastData);
        lenOfArr = dim(1,2);

        for i = 1:lenOfArr
            if lastData(1,i) == 1
                ss_copies_count  = ss_copies_count  + lastData(3,i);
            elseif lastData(1,i) == 2
                pds_copies_count = pds_copies_count + lastData(3,i); 
            elseif lastData(1,i) == 3
                sl_copies_count  = sl_copies_count  + lastData(3,i);
            else
                t_copies_count   = t_copies_count   + lastData(3,i);
            end
        end

        lastAmpliconsLen = [lastAmpliconsLen,lastData(2:3,:)];
        lengths          = lastAmpliconsLen';
        lengths          = sortrows(lengths,1);
        [r,c]            = size(lengths);

        uniqueLens       = unique(lengths(:,1));
        countOfLen       = hist(lengths(:,1),uniqueLens);

        len = 1;
        m   = 1;

        while true
            lengthData(m,1) = lengths(len,1);
            lengthData(m,2) = sum(lengths(len:(len + countOfLen(1,m) - 1),2));
            len = len + countOfLen(1,m);
            m  = m + 1; 
            if len > r
                break;
            end
        end

        lastAmpliconsLen = lengthData';
        
    end
    
    % STEP 4 : EXTRACTING AND ADDING THE DATA INTO SLPACKAGEBUFFER1 AND PDSPACKAGEBUFFER
    % 'y' for loop deals with first the amplicons from pathway1 and then in the second iteration with the amplicons from pathway2
    for y = 1:numAmpliconInput
        
        % 4.1 Sending the SL data from PDS1 and PDS2 into SL package buffer
        if y == 1
            
            dataSLfPDS1  = dataFun{11};
            dim          = size(dataSLfPDS1);
            numRows      = dim(1,1);
            numCols      = dim(1,2);
            numIter      = numCols/2;
            
            shiftCol = 0;

            for z = 1:numIter
                slPackage{1} = dataFun{9};
                slTimeArr    = dataSLfPDS1(1:numRows,1 + shiftCol);
                slPackage{2} = slTimeArr;
                slPackage{3} = dataSLfPDS1(1:numRows,2 + shiftCol);
                if slTimeArr(1,1) <= simTime
                    slPackageBuffer2{P1bufferSL2,1} = slPackage;
                    P1bufferSL2  = P1bufferSL2 + 1;
                end    
                shiftCol = shiftCol + 2;
            end
            
            dataSLfPDS2  = dataFun{12};
            dim          = size(dataSLfPDS2);
            numRows      = dim(1,1);
            numCols      = dim(1,2); 
            numIter      = numCols/2;
            
            shiftCol = 0;

            for z = 1:numIter
                slPackage{1} = dataFun{10};
                slTimeArr    = dataSLfPDS2(1:numRows,1 + shiftCol);
                slPackage{2} = slTimeArr;
                slPackage{3} = dataSLfPDS2(1:numRows,2 + shiftCol);
                if slTimeArr(1,1) <= simTime
                    slPackageBuffer2{P1bufferSL2,1} = slPackage;
                    P1bufferSL2  = P1bufferSL2 + 1;
                end
                shiftCol = shiftCol + 2;
            end
                
        else
            
            dataSLfPDS1  = dataFun{26};
            dim          = size(dataSLfPDS1);
            numRows      = dim(1,1);
            numCols      = dim(1,2);
            numIter      = numCols/2;
            
            shiftCol = 0;
             
            for z = 1:numIter              
                slPackage{1} = dataFun{24};
                slTimeArr    = dataSLfPDS1(1:numRows,1 + shiftCol);
                slPackage{2} = slTimeArr;
                slPackage{3} = dataSLfPDS1(1:numRows,2 + shiftCol);
                if slTimeArr(1,1) <= simTime
                    slPackageBuffer2{P2bufferSL2,2} = slPackage;
                    P2bufferSL2  = P2bufferSL2 + 1;
                end
                shiftCol = shiftCol + 2;
            end
            
            dataSLfPDS2  = dataFun{27};
            dim          = size(dataSLfPDS2);
            numRows      = dim(1,1);
            numCols      = dim(1,2);
            numIter      = numCols/2;
            
            shiftCol = 0;

            for z = 1:numIter
                slPackage{1} = dataFun{25};
                slTimeArr    = dataSLfPDS2(1:numRows,1 + shiftCol);
                slPackage{2} = slTimeArr;
                slPackage{3} = dataSLfPDS2(1:numRows,2 + shiftCol);
                if slTimeArr(1,1) <= simTime
                    slPackageBuffer2{P2bufferSL2,2} = slPackage;
                    P2bufferSL2  = P2bufferSL2 + 1;
                end
                shiftCol = shiftCol + 2;               
            end
            
        end
        
        if y == 1
            childPDSnum = dataFun{13};
        else
            childPDSnum = dataFun{28};
        end
        
        % 4.2 Sending the PDS data from SS1 and SS2 into PDS package buffer
        % z -> Deals with the PDS from SS1 and SS2
        % w -> PDS from SS1 and SS2 are produced at different times according to the starting initial SS time
        % v -> takes care of how many child PDS are produced from SS1 and SS2
        for z = 1:2
            
            if y == 1
                if z == 1
                    strPDSfSS     = dataFun{1};
                    dataPDSfSS    = dataFun{2};
                    numOfChildPDS = childPDSnum(1,1);
                else
                    strPDSfSS     = dataFun{5};
                    dataPDSfSS    = dataFun{6};
                    numOfChildPDS = childPDSnum(1,3);
                end
            else
                if z == 1
                    strPDSfSS      = dataFun{16};
                    dataPDSfSS     = dataFun{17};
                    numOfChildPDS  = childPDSnum(1,1);
                else
                    strPDSfSS      = dataFun{20};
                    dataPDSfSS     = dataFun{21};
                    numOfChildPDS  = childPDSnum(1,3); 
                end
            end
            
            % If its empty that means no PDS was formed from SS
            % If its not empty then, childPDSnum will tell you how many child PDS were formed from that SS
            if (isempty(strPDSfSS) == 0)      
                dim       = size(dataPDSfSS);
                numOfCols = dim(1,2);
                numOfIter = numOfCols/(numOfChildPDS + 1);
                shiftCol  = 0;    

                % For loop 'w' deals with the PDS amplicons generated at different times from different initial SS times 
                for w = 1:numOfIter
                    % For loop 'v' deals with the different child PDS amplicons
                    for v = 1:numOfChildPDS
                        if y == 1
                            pdsPackage1{1} = strPDSfSS{1,v};
                            pdsPackage1{2} = strPDSfSS{2,v};
                            pdsTimeArr     = dataPDSfSS(1:end,(v + shiftCol));
                            pdsPackage1{3} = pdsTimeArr;
                            pdsPackage1{4} = dataPDSfSS(1:end,(numOfChildPDS + 1 + shiftCol));
                            if pdsTimeArr(1,1) <= simTime
                                pdsPackageBuffer{P1bufferPDS,1} = pdsPackage1;
                                P1bufferPDS = P1bufferPDS + 1;
                            end
                        else
                            pdsPackage2{1} = strPDSfSS{1,v};
                            pdsPackage2{2} = strPDSfSS{2,v};
                            pdsTimeArr     = dataPDSfSS(1:end,(v + shiftCol));
                            pdsPackage2{3} = pdsTimeArr;  
                            pdsPackage2{4} = dataPDSfSS(1:end,(numOfChildPDS + 1 + shiftCol));
                            if pdsTimeArr(1,1) <= simTime
                                pdsPackageBuffer{P2bufferPDS,2} = pdsPackage2; 
                                P2bufferPDS = P2bufferPDS + 1;
                            end
                        end
                    end
                    shiftCol  = shiftCol  + numOfChildPDS + 1;
                end
            end
            
        end
        
        % 4.3 Sending the PDS data from PDS1 and PDS2 into PDS package buffer
        for z = 1:2
            
            if y == 1
                if z == 1
                    strPDSfPDS    = dataFun{3};
                    dataPDSfPDS   = dataFun{4};
                    numOfChildPDS = childPDSnum(1,2);
                else
                    strPDSfPDS    = dataFun{7};
                    dataPDSfPDS   = dataFun{8};
                    numOfChildPDS = childPDSnum(1,4);
                end
            else
                if z == 1
                    strPDSfPDS    = dataFun{18};
                    dataPDSfPDS   = dataFun{19};
                    numOfChildPDS = childPDSnum(1,2);
                else
                    strPDSfPDS    = dataFun{22};
                    dataPDSfPDS   = dataFun{23};
                    numOfChildPDS = childPDSnum(1,4); 
                end
            end
            
            % This checks whether a PDS was formed from a PDS
            if (isempty(strPDSfPDS) == 0)
                
                dim       = size(dataPDSfPDS);
                numOfCols = dim(1,2);
                % The + 1 takes care of the copies column also
                numOfIter = numOfCols/(numOfChildPDS + 1);
                shiftCol  = 0;

                % For loop 'w' deals with the PDS amplicons generated at different times from different initial SS times from parent PDS belonging to the highway
                for w = 1:numOfIter
                    % For loop 'v' deals with the different child PDS amplicons
                    for v = 1:numOfChildPDS
                        if y == 1
                            pdsPackage1{1} = strPDSfPDS{1,v};
                            pdsPackage1{2} = strPDSfPDS{2,v};
                            pdsTimeArr     = dataPDSfPDS(1:end,(v + shiftCol));
                            pdsPackage1{3} = pdsTimeArr;
                            pdsPackage1{4} = dataPDSfPDS(1:end,(numOfChildPDS + 1 + shiftCol));
                            if pdsTimeArr(1,1) <= simTime
                                pdsPackageBuffer{P1bufferPDS,1} = pdsPackage1;
                                P1bufferPDS = P1bufferPDS + 1;
                            end
                        else
                            pdsPackage2{1} = strPDSfPDS{1,v};
                            pdsPackage2{2} = strPDSfPDS{2,v};
                            pdsTimeArr     = dataPDSfPDS(1:end,(v + shiftCol));
                            pdsPackage2{3} = pdsTimeArr;
                            pdsPackage2{4} = dataPDSfPDS(1:end,(numOfChildPDS + 1 + shiftCol));
                            if pdsTimeArr(1,1) <= simTime
                                pdsPackageBuffer{P2bufferPDS,2} = pdsPackage2;
                                P2bufferPDS = P2bufferPDS + 1;
                            end
                        end
                    end
                    shiftCol  = shiftCol  + numOfChildPDS + 1;
                end
            end
            
        end
            
    % --- END OF FOR LOOP ---     
    end
 
% --- END OF LOOP CONSIDERING DATA FROM 8 SS AMPLICONS ---     
end

end
%----------------------------------- END OF SECTION 1 -----------------------------------



% --------------------------------------- SECTION 2 --------------------------------------

% This gets activated when space in the pdsPackageBuffer is running low
if increaseNumFun == 1
    numPackages = 40;
else
    numPackages = 8;
end
dim = size(pdsPackageBuffer);
if dim(1,1) <= numPackages
    transferCell     = pdsPackageBuffer(1:end,:); 
    pdsPackageBuffer = cell(replenishRowAmt,2);
    pdsPackageBuffer(1:dim(1,1),:) = transferCell;
end

% Finding how many packages in the PDS buffer
findEmptySlot   = find(cellfun(@isempty,pdsPackageBuffer),1);
lastRowPDS      = findEmptySlot - 1;
if isempty(lastRowPDS) == 1
    lastRowPDS = dim(1,1);
end

if lastRowPDS ~= 0

if increaseNumFun == 1
    
    % Running 40 functions simultaneuosly
    funcs     = {@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS};
    arguments = {pdsPackageBuffer{1,1},pdsPackageBuffer{1,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{2,1},pdsPackageBuffer{2,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{3,1},pdsPackageBuffer{3,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{4,1},pdsPackageBuffer{4,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{5,1},pdsPackageBuffer{5,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{6,1},pdsPackageBuffer{6,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{7,1},pdsPackageBuffer{7,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{8,1},pdsPackageBuffer{8,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen
                 pdsPackageBuffer{9,1},pdsPackageBuffer{9,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{10,1},pdsPackageBuffer{10,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{11,1},pdsPackageBuffer{11,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{12,1},pdsPackageBuffer{12,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{13,1},pdsPackageBuffer{13,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{14,1},pdsPackageBuffer{14,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{15,1},pdsPackageBuffer{15,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{16,1},pdsPackageBuffer{16,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{17,1},pdsPackageBuffer{17,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{18,1},pdsPackageBuffer{18,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{19,1},pdsPackageBuffer{19,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{20,1},pdsPackageBuffer{20,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{21,1},pdsPackageBuffer{21,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{22,1},pdsPackageBuffer{22,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{23,1},pdsPackageBuffer{23,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{24,1},pdsPackageBuffer{24,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{25,1},pdsPackageBuffer{25,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{26,1},pdsPackageBuffer{26,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{27,1},pdsPackageBuffer{27,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{28,1},pdsPackageBuffer{28,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{29,1},pdsPackageBuffer{29,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{30,1},pdsPackageBuffer{30,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{31,1},pdsPackageBuffer{31,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{32,1},pdsPackageBuffer{32,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{33,1},pdsPackageBuffer{33,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{34,1},pdsPackageBuffer{34,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{35,1},pdsPackageBuffer{35,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{36,1},pdsPackageBuffer{36,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{37,1},pdsPackageBuffer{37,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{38,1},pdsPackageBuffer{38,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{39,1},pdsPackageBuffer{39,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{40,1},pdsPackageBuffer{40,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen};
    solutions = cell(40,1); 
    
    % Clearing the entires that have been used up
    pdsPackageBuffer(1:40,:) = [];
    P1bufferPDS = P1bufferPDS - 40;
    P2bufferPDS = P2bufferPDS - 40;
    if P1bufferPDS <= 0
        P1bufferPDS = 1;
    end
    if P2bufferPDS <= 0
        P2bufferPDS = 1;
    end
    
    % Here the 28 functions run parallely
    parfor ii = 1:40
        solutions{ii} = funcs{ii}(arguments{ii,:});
    end
    
    numOfIterX = 40;
    
else
    
    funcs     = {@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS,@childPDS};
    arguments = {pdsPackageBuffer{1,1},pdsPackageBuffer{1,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{2,1},pdsPackageBuffer{2,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{3,1},pdsPackageBuffer{3,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{4,1},pdsPackageBuffer{4,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{5,1},pdsPackageBuffer{5,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{6,1},pdsPackageBuffer{6,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{7,1},pdsPackageBuffer{7,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 pdsPackageBuffer{8,1},pdsPackageBuffer{8,2},dataMatBlank,simTime,timeInt,rate,probSL,patchLen};
    solutions = cell(8,1);
    
    % Deleting the packages that have been sent into the parfor
    pdsPackageBuffer(1:8,:) = [];
    P1bufferPDS = P1bufferPDS - 8;
    P2bufferPDS = P2bufferPDS - 8;
    if P1bufferPDS <= 0
        P1bufferPDS = 1;
    end
    if P2bufferPDS <= 0
        P2bufferPDS = 1;
    end
    
    % Here the 8 functions run parallely
    parfor ii = 1:8
        solutions{ii} = funcs{ii}(arguments{ii,:});
    end
    
    numOfIterX = 8;
    
end

for x = 1:numOfIterX
    
    dataFun = solutions{x};
    
    % STEP 1 : Determining the length of dataFun
    if length(dataFun) == 6
        numAmpliconInput = 1;
        lastData         = dataFun{5};
        addData          = dataFun{6};
    elseif length(dataFun) == 12
        numAmpliconInput = 2;
        lastData         = dataFun{11};
        addData          = dataFun{12};
    else
        continue;
    end
    
    % STEP 2 : Adding the data into dataMat
    for y = 2:6
        dataMat(1:end,y) = dataMat(1:end,y) + addData(1:end,y);
    end
    
    % STEP 3 : Adding the amplicons formed at the end into the lastAmpliconsLen
    if isempty(lastData) == 0
        
        dim = size(lastData);
        lenOfArr = dim(1,2);

        for i = 1:lenOfArr
            if lastData(1,i) == 1
                ss_copies_count  = ss_copies_count  + lastData(3,i);
            elseif lastData(1,i) == 2
                pds_copies_count = pds_copies_count + lastData(3,i); 
            elseif lastData(1,i) == 3
                sl_copies_count  = sl_copies_count  + lastData(3,i);
            else
                t_copies_count   = t_copies_count   + lastData(3,i);
            end
        end
        
        lastAmpliconsLen = [lastAmpliconsLen,lastData(2:3,:)]; %#ok<*AGROW>
        lengths          = lastAmpliconsLen';
        lengths          = sortrows(lengths,1);
        [r,c]            = size(lengths);

        uniqueLens       = unique(lengths(:,1));
        countOfLen       = hist(lengths(:,1),uniqueLens);

        len = 1;
        m   = 1;

        while true
            lengthData(m,1) = lengths(len,1);
            lengthData(m,2) = sum(lengths(len:(len + countOfLen(1,m) - 1),2));
            len = len + countOfLen(1,m);
            m  = m + 1; 
            if len > r
                break;
            end
        end

        lastAmpliconsLen = lengthData';
    end
    
    % For loop 'y' deals with the amplicons from pathway 1 and pathway 2
    for y = 1:numAmpliconInput
    
        % STEP 4 : Sending the SS data into sspackageBuffer
        if y == 1
            strSSfPDS  = dataFun{1};
            dataSSfPDS = dataFun{2};

            ssPackage{1} = strSSfPDS{1,1};
            ssPackage{2} = strSSfPDS{2,1};
            ssTimeArr    = dataSSfPDS(1:end,1);
            ssPackage{3} = ssTimeArr;
            ssPackage{4} = dataSSfPDS(1:end,2);
            if ssTimeArr(1,1) <= simTime
                ssPackageBuffer{P1bufferSS,1} = ssPackage;
                P1bufferSS = P1bufferSS + 1;
            end
        else
            strSSfPDS  = dataFun{7};
            dataSSfPDS = dataFun{8}; 

            ssPackage{1} = strSSfPDS{1,1};
            ssPackage{2} = strSSfPDS{2,1};
            ssTimeArr     = dataSSfPDS(1:end,1);
            ssPackage{3} = ssTimeArr;
            ssPackage{4} = dataSSfPDS(1:end,2);
            if ssTimeArr(1,1) <= simTime
                ssPackageBuffer{P2bufferSS,2} = ssPackage;
                P2bufferSS = P2bufferSS + 1;
            end
        end
    
        % STEP 5 : Sending the SL data into slpackageBuffer2
        if y == 1
            dsSLfPDS   = dataFun{3};
            dataSLfPDS = dataFun{4};

            slPackage{1} = dsSLfPDS;
            slTimeArr    = dataSLfPDS(1:end,1);
            slPackage{2} = slTimeArr;
            slPackage{3} = dataSLfPDS(1:end,2);
            if slTimeArr(1,1) <= simTime
                slPackageBuffer1{P1bufferSL1,1} = slPackage;
                P1bufferSL1 = P1bufferSL1 + 1;
            end
        else
            dsSLfPDS   = dataFun{9};
            dataSLfPDS = dataFun{10};

            slPackage{1} = dsSLfPDS;
            slTimeArr    = dataSLfPDS(1:end,1);
            slPackage{2} = slTimeArr;
            slPackage{3} = dataSLfPDS(1:end,2);
            if slTimeArr(1,1) <= simTime
                slPackageBuffer1{P2bufferSL1,2} = slPackage;
                P2bufferSL1 = P2bufferSL1 + 1;
            end
        end
    
    end
    
end
    
end
% ----------------------------------- END OF SECTION 2 -----------------------------------



% -------------------------------------- SECTION 3 ---------------------------------------

% Here we deal with two different types of slPackageBuffers. 
% 1) slPackageBuffer1 gets its data packages from the childPDS function
% 2) slPackageBuffer2 gets its data packages from the ssLoop function

for m = 1:2

    if increaseNumFun == 1
        % Initializing the storage containers to send the data to parfor loop
        dataContainer1 = cell(40,1);
        dataContainer2 = cell(40,1);
    else
        % Initializing the storage containers to send the data to parfor loop
        dataContainer1 = cell(8,1);
        dataContainer2 = cell(8,1);
    end

%  ----- THIS IS FOR DATA OF SL PRODUCED FROM THE childPDS -----
if m == 1
    
    if increaseNumFun == 1
        numOfIterY = 40;
    else
        numOfIterY = 8;
    end
 
    dim = size(slPackageBuffer1);
    if dim(1,1) <= numOfIterY
        transferCell     = slPackageBuffer1(1:end,:); 
        slPackageBuffer1 = cell(replenishRowAmt,2);
        slPackageBuffer1(1:dim(1,1),:) = transferCell;
    end
    
    % Finding out the number of packages left to analyze in slPackageBuffer1
    findEmptySlot1  = find(cellfun(@isempty,slPackageBuffer1),1);
    lastRow1        = findEmptySlot1 - 1;
    if isempty(lastRow1) == 1
        lastRow1 = dim(1,1);
    end
    
    % Here we check that if slPackageBuffer1 is empty then continue with the next iteration
    if lastRow1 == 0
        continue;
    end

    for y = 1:numOfIterY  
        % Sending the first 8/28 data packages into the data containers 
        dataContainer1{y,1} = slPackageBuffer1{y,1};        
        dataContainer2{y,1} = slPackageBuffer1{y,2};
    end
    
    if numOfIterY == 40
        slPackageBuffer1(1:40,:) = [];
        P1bufferSL1 = P1bufferSL1 - 40;
        P2bufferSL1 = P2bufferSL1 - 40;
    else
        slPackageBuffer1(1:8,:) = [];
        P1bufferSL1 = P1bufferSL1 - 8;
        P2bufferSL1 = P2bufferSL1 - 8;
    end
    
    if P1bufferSL1 <= 0
        P1bufferSL1 = 1;
    end
    if P2bufferSL1 <= 0
        P2bufferSL1 = 1;
    end   

%  ----- THIS IS FOR DATA OF SL PRODUCED FROM THE ssLOOP -----
else
     
    if increaseNumFun == 1
        numOfIterY = 40;
    else
        numOfIterY = 8;
    end
    
    dim = size(slPackageBuffer2);
    if dim(1,1) <= numOfIterY
        transferCell     = slPackageBuffer2(1:end,:); 
        slPackageBuffer2 = cell(replenishRowAmt,2);
        slPackageBuffer2(1:dim(1,1),:) = transferCell;
    end
    
    % Finding out the number of packages left to analyze in slPackageBuffer2
    findEmptySlot2  = find(cellfun(@isempty,slPackageBuffer2),1);
    lastRow2        = findEmptySlot2 - 1;
    if isempty(lastRow2) == 1
        lastRow2 = dim(1,1);
    end
    
    % Here we check that if slPackageBuffer2 is empty then continue with the next iteration
    if lastRow2 == 0
        break;
    end
    
    for y = 1:numOfIterY  
        % Sending the first 8/28 data packages into the data containers 
        dataContainer1{y,1} = slPackageBuffer2{y,1};  
        dataContainer2{y,1} = slPackageBuffer2{y,2};
    end
    
    if numOfIterY == 40
        slPackageBuffer2(1:40,:) = [];
        P1bufferSL2 = P1bufferSL2 - 40;
        P2bufferSL2 = P2bufferSL2 - 40;
    else
        slPackageBuffer2(1:8,:) = [];
        P1bufferSL2 = P1bufferSL2 - 8;
        P2bufferSL2 = P2bufferSL2 - 8;
    end
    
    if P1bufferSL2 <= 0
        P1bufferSL2 = 1;
    end
    if P2bufferSL2 <= 0
        P2bufferSL2 = 1;
    end
    
% --- END OF FOR LOOP 'X' ----
end

% Sending the data containers filled with the SL data into the SL functions

if increaseNumFun == 1
    
    funcs     = {@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy};
    arguments = {dataContainer1{1,1},dataContainer2{1,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{2,1},dataContainer2{2,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{3,1},dataContainer2{3,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{4,1},dataContainer2{4,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{5,1},dataContainer2{5,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{6,1},dataContainer2{6,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{7,1},dataContainer2{7,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{8,1},dataContainer2{8,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{9,1},dataContainer2{9,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{10,1},dataContainer2{10,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{11,1},dataContainer2{11,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{12,1},dataContainer2{12,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{13,1},dataContainer2{13,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{14,1},dataContainer2{14,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{15,1},dataContainer2{15,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{16,1},dataContainer2{16,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{17,1},dataContainer2{17,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{18,1},dataContainer2{18,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{19,1},dataContainer2{19,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{20,1},dataContainer2{20,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{21,1},dataContainer2{21,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{22,1},dataContainer2{22,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{23,1},dataContainer2{23,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{24,1},dataContainer2{24,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{25,1},dataContainer2{25,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{26,1},dataContainer2{26,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{27,1},dataContainer2{27,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{28,1},dataContainer2{28,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{29,1},dataContainer2{29,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{30,1},dataContainer2{30,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{31,1},dataContainer2{31,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{32,1},dataContainer2{32,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{33,1},dataContainer2{33,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{34,1},dataContainer2{34,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{35,1},dataContainer2{35,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{36,1},dataContainer2{36,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{37,1},dataContainer2{37,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{38,1},dataContainer2{38,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{39,1},dataContainer2{39,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{40,1},dataContainer2{40,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen};
    solutions = cell(40,1);

    parfor ii = 1:40
        solutions{ii} = funcs{ii}(arguments{ii,:});
    end
    
    numOfIterX = 40;
    
else
    
    funcs     = {@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy,@slHwy};
    arguments = {dataContainer1{1,1},dataContainer2{1,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{2,1},dataContainer2{2,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{3,1},dataContainer2{3,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{4,1},dataContainer2{4,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{5,1},dataContainer2{5,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{6,1},dataContainer2{6,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{7,1},dataContainer2{7,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen;
                 dataContainer1{8,1},dataContainer2{8,1},dataMatBlank,simTime,timeInt,rate,probSL,patchLen};
    solutions = cell(8,1);

    parfor ii = 1:8
        solutions{ii} = funcs{ii}(arguments{ii,:});
    end
    
    numOfIterX = 8;
    
end

% Extracting the data from solutions cell
for  x = 1:numOfIterX
    
    dataFun = solutions{x};
   
    % STEP 1 : DETERMINING THE LENGTH OF DATAFUN
    % This is done to find out how many amplicons were input(0,1,2). This helps to figure out how many times for loop 'x' is to be cycled
    % If length is 7, then one amplicon was input
    if length(dataFun) == 7 
        numAmpliconInput = 1;
        lastData         = dataFun{6};
        addData          = dataFun{7};
    % If length is 12, then two amplicons were input
    elseif length(dataFun) == 14
        numAmpliconInput = 2;
        lastData         = dataFun{13};
        addData          = dataFun{14};
    % This corresponds to no amplicons were input. Then proceed with the next iteration
    else
        continue;
    end
    
    % NOTE : In a function addData and lastData are common for both sister amplicons
    
    % STEP 2 : ADDING THE addData FROM FUNCTION TO THE COMMON dataMat
    for y = 2:6
        dataMat(1:end,y) = dataMat(1:end,y) + addData(1:end,y);
    end
    
    % STEP 3 : ADDING THE lastData FROM FUNCTION TO THE COMMON lastAmpliconsData
    if isempty(lastData) == 0
        
        dim = size(lastData);
        lenOfArr = dim(1,2);

        for i = 1:lenOfArr
            if lastData(1,i) == 1
                ss_copies_count  = ss_copies_count  + lastData(3,i);
            elseif lastData(1,i) == 2
                pds_copies_count = pds_copies_count + lastData(3,i); 
            elseif lastData(1,i) == 3
                sl_copies_count  = sl_copies_count  + lastData(3,i);
            else
                t_copies_count   = t_copies_count   + lastData(3,i);
            end
        end

        lastAmpliconsLen = [lastAmpliconsLen,lastData(2:3,:)];
        lengths          = lastAmpliconsLen';
        lengths          = sortrows(lengths,1);
        [r,c]            = size(lengths);

        uniqueLens       = unique(lengths(:,1));
        countOfLen       = hist(lengths(:,1),uniqueLens);

        len = 1;
        k   = 1;

        while true
            lengthData(k,1) = lengths(len,1);
            lengthData(k,2) = sum(lengths(len:(len + countOfLen(1,k) - 1),2));
            len = len + countOfLen(1,k);
            k   = k + 1; 
            if len > r
                break;
            end
        end
        lastAmpliconsLen = lengthData';

    end
    
    % STEP 4 : EXTRACTING THE SS AND PDS DATA 
    for y = 1:numAmpliconInput
         
        % 4.1 UNPACKING THE SS DATA
        if y == 1
            strSSfPDS   = dataFun{1};
            dataSSfPDS  = dataFun{2};
        else
            strSSfPDS   = dataFun{8};
            dataSSfPDS  = dataFun{9};
        end

        if isempty(strSSfPDS) == 0
            dim         = size(dataSSfPDS);
            numOfCol    = dim(1,2);
            numOfIter   = numOfCol/2;
            shiftCol    = 0;

            % Creating the SS buffer array
            for z = 1:numOfIter
                if y == 1
                    ssPackage1{1}        = strSSfPDS{1,z};
                    ssPackage1{2}        = strSSfPDS{2,z};
                    ssTimeArr            = dataSSfPDS(1:end,1 + shiftCol);
                    ssPackage1{3}        = ssTimeArr;
                    ssPackage1{4}        = dataSSfPDS(1:end,2 + shiftCol);
                    if ssTimeArr(1,1) <= simTime
                        ssPackageBuffer{P1bufferSS,1} = ssPackage1; 
                        P1bufferSS = P1bufferSS + 1;
                    end
                else
                    ssPackage2{1}        = strSSfPDS{1,z};
                    ssPackage2{2}        = strSSfPDS{2,z};
                    ssTimeArr            = dataSSfPDS(1:end,1 + shiftCol);
                    ssPackage2{3}        = ssTimeArr;
                    ssPackage2{4}        = dataSSfPDS(1:end,2 + shiftCol);
                    if ssTimeArr(1,1) <= simTime
                        ssPackageBuffer{P2bufferSS,2} = ssPackage2;
                        P2bufferSS = P2bufferSS + 1;
                    end
                end
                shiftCol = shiftCol + 2;
            end
        end

        % 4.2 UNPACKING THE PDS DATA
        if y == 1
            strPDSfPDS  = dataFun{3};
            dataPDSfPDS = dataFun{4};
            childPDSnum = dataFun{5};
        else
            strPDSfPDS  = dataFun{10};
            dataPDSfPDS = dataFun{11};
            childPDSnum = dataFun{12};
        end
        
        % The childPDSnum tells how many child PDS were produced from the parent PDS
        % length of this array tells how many iterations muct be done in for loop 'z'
        lengthOfArr = length(childPDSnum);
        shiftCol1   = 0;
        shiftCol2   = 0;

        if isempty(strPDSfPDS) == 0
            for z = 1:lengthOfArr
                numOfChildPDS   = childPDSnum(1,z);
                pdsCopiesArr    = dataPDSfPDS(1:end,(numOfChildPDS + 1 + shiftCol2));
                for w = 1:numOfChildPDS
                    if y == 1
                        pdsPackage1{1}          = strPDSfPDS{1,w + shiftCol1};
                        pdsPackage1{2}          = strPDSfPDS{2,w + shiftCol1};
                        pdsTimeArr              = dataPDSfPDS(1:end,w + shiftCol2);
                        pdsPackage1{3}          = pdsTimeArr;
                        pdsPackage1{4}          = pdsCopiesArr;
                        if pdsTimeArr(1,1) <= simTime
                            pdsPackageBuffer{P1bufferPDS,1} = pdsPackage1;
                            P1bufferPDS = P1bufferPDS + 1;
                        end
                    else
                        pdsPackage2{1}          = strPDSfPDS{1,w + shiftCol1};
                        pdsPackage2{2}          = strPDSfPDS{2,w + shiftCol1};
                        pdsTimeArr              = dataPDSfPDS(1:end,w + shiftCol2);
                        pdsPackage2{3}          = pdsTimeArr;
                        pdsPackage2{4}          = pdsCopiesArr;
                        if pdsTimeArr(1,1) <= simTime
                            pdsPackageBuffer{P2bufferPDS,2} = pdsPackage2;
                            P2bufferPDS = P2bufferPDS + 1;
                        end
                    end
                end
                shiftCol1 = shiftCol1 + numOfChildPDS;
                shiftCol2 = shiftCol2 + numOfChildPDS + 1;
            end
        end
        
    end    
% --- END OF FOR LOOP 'x' ---  
end

end

% EXITING FROM THE WHILE LOOP
if (lastRowSS == 0) && (lastRowPDS == 0) && (lastRow1 == 0) && (lastRow2 == 0)
    break;
end

% --- END OF WHILE (GOES THROUGH ALL 3 SECTIONS) ---
end
% ----------------------------------- END OF SECTION 3 -----------------------------------

% --------------------------------- DISPLAYING THE DATA ----------------------------------

disp("NUMBER OF COPIES OF SS AT THE END OF SIMTIME  : ");
disp(ss_copies_count);
disp("NUMBER OF COPIES OF PDS AT THE END OF SIMTIME : ");
disp(pds_copies_count);
disp("NUMBER OF COPIES OF SL AT THE END OF SIMTIME  : ");
disp(sl_copies_count);
disp("NUMBER OF COPIES OF T AT THE END OF SIMTIME   : ");
disp(t_copies_count);

totalCopies = sum(lastAmpliconsLen(2,:));

disp("TOTAL NUMBER OF COPIES                        : ");
disp(totalCopies);

% Finally we add the cumulative nucletoides to see how the number of nucleotides vary with time
speciesPlot(1,1) = 0;
speciesPlot(1,2) = copiesOfDNA;
for i = 2:rowsDataMat
    speciesPlot(i,1) = dataMat(i,1);
    speciesPlot(i,2) = dataMat(i,2) + dataMat(i,3) + dataMat(i,4) + dataMat(i,5);
    dataMat(i,6)     = dataMat(i - 1,6) + dataMat(i,6);
end

xAxis_time   = speciesPlot(1:rowsDataMat,1);
yAxis_copies = speciesPlot(1:rowsDataMat,2);
yAxis_nucleo = dataMat(1:rowsDataMat,6);

save(specifyFileName);

% % DISPLAYING THE SIMULATED GEL
% 
% lastAmpliconsLen = lastAmpliconsLen';
% largestLength    = lastAmpliconsLen(end,1);
% lastAmpliconsLen = sortrows(lastAmpliconsLen,-2);
% largestNumber    = lastAmpliconsLen(1,2);
% 
% %y axis maximum value
% ymax = largestLength + 50;
% %x axis maximum value
% xmax = 15;
% %Defining the plot area
% axis    = [0 xmax 0 ymax];
% mainRec = rectangle('Position',[0 0 xmax ymax]);
% mainRec.FaceColor = [0 0 0];
% 
% %width and height of rectangle
% width  = xmax;
% height = 10; 
% 
% % Size of the lastAmpliconsLen matrix
% dim = size(lastAmpliconsLen);
% r   = dim(1,1);
% 
% for i = 1:r
%     length = lastAmpliconsLen(i,1);
%     fc     = lastAmpliconsLen(i,2)/largestNumber;
%     rectangle('Position',[0 length width height],'FaceColor',[fc fc fc])
% end





