function [dataFun] = ssLoop(ssData1,ssData2,dataMat,simTime,timeInt,rate,probSL,patchLen)

% Outline of this loop  is 
% SS1 -> SL -> PDS1 -> SS2 -> SL -> PDS2 -> SS1
% A unique fact about this function : 
% 1. In ssTimeArr we consider each inidividual time as SS1 is being replenished after a fixed period of time. Hence the time arrays for
%    formation of SL and PDS would be different for different starting values of ssformationTime.
% 2. All the amplicons participating in the loop have the same length.

% Outer for loop takes care of the sister amplicons formed from pathways 1 and 2 steming from the dumbell loop
% EX : [B2 F2 B2C] and [F2 B2 F2C]

% CASE 1 : NO AMPLICON IS INPUT
if (isempty(ssData1) && isempty(ssData2))
    dataFun = cell(0,0);
    return 
% CASE 2 : AMPLICON FROM PATHWAY 1 IS INPUT
elseif (isempty(ssData1))
    ss           = ssData2{1};
    css          = ssData2{2};
    ssTimeArr    = ssData2{3};
    ssCopiesArr  = ssData2{4};
    breakTheLoop = 1;
% CASE 3 : AMPLICON FROM PATHWAY 2 IS INPUT
elseif (isempty(ssData2))
    ss           = ssData1{1};
    css          = ssData1{2};
    ssTimeArr    = ssData1{3};
    ssCopiesArr  = ssData1{4};
    breakTheLoop = 1;
% CASE 4 : BOTH AMPLICONS INPUT
else
    breakTheLoop = 0;
end

% This stores the length of the last Amplicons
lastAmpliconsLen = zeros(0,0);
la = 1;

for x = 1:2
    
if breakTheLoop == 0
    if x == 1
        ss           = ssData1{1};
        css          = ssData1{2};
        ssTimeArr    = ssData1{3};
        ssCopiesArr  = ssData1{4};
    else
        ss           = ssData2{1};
        css          = ssData2{2};
        ssTimeArr    = ssData2{3};
        ssCopiesArr  = ssData2{4};
    end
end

% Finding the length of the ss array
% Here by length we mean the number of B and F patches in SS
lenOfArrSS = length(ss);

if rem(lenOfArrSS,2) == 0
    bcount = lenOfArrSS/2;
    fcount = lenOfArrSS/2;
    
    % For the participants in the loop, only SS and SL have the same length
    lengthSS  = [bcount,fcount,(bcount + fcount - 1)]*patchLen;% This gives the number of nucleotides present in SS
    lengthSL  = lengthSS;
    naPDS     = lengthSS;% PDS formed from SL
else
    bcount = ceil(lenOfArrSS/2);
    fcount = lenOfArrSS - bcount;
    if ss(1,end) == "B2c"
        lengthSS = [bcount,fcount,(bcount + fcount - 1)]*patchLen;
    else
        lengthSS = [fcount,bcount,(bcount + fcount - 1)]*patchLen;
    end
    lengthSL  = lengthSS;
    naPDS     = lengthSS;
end

% Creating the cells and arrays that would be used in this function
% The below containers contain the child PDS produced from SS/PDS
%     -- PDSfSS1  --
strPDSfSS1    = cell(0,0);
dataPDSfSS1   = zeros(0,0);
%     -- PDSfPDS1 --
strPDSfPDS1   = cell(0,0);
dataPDSfPDS1  = zeros(0,0);
%     -- PDSfSS2  --
strPDSfSS2    = cell(0,0);
dataPDSfSS2   = zeros(0,0);
%     -- PDSfPDS2 --
strPDSfPDS2   = cell(0,0);
dataPDSfPDS2  = zeros(0,0);
%     -- SLfPDS1  --
dsSLfPDS1     = strings(0,0);
dsSLfPDS2     = strings(0,0);
dataSLfPDS1   = zeros(0,0);
dataSLfPDS2   = zeros(0,0);

% The below array will contain the number of child PDS that can be produced from a parent PDS  
% (provided the number of copies of child PDS produced are greater than 0)
childPDSnum   = zeros(1,4);

% Inner for loop splits the ss loop into two parts and takes care of them in separate iterations
for y = 1:2
% Consider SS1 -> SL1 -> PDS1 for iteration 1
% SS2 -> SL2 -> PDS2 for iteration 2 (Here SS2 is just the complementary version of SS1) 

% ------------------------------------------------------------------------

% SS1/SS2
if y == 1
    extensionTimeArr(1,1) = 0;
    numOfCopiesSS         = ssCopiesArr(1,1);  
else
% STEP 1 : SENDING SS2 INFO INTO dataMat
    prevss = ss;
    ss     = css;
    css    = prevss;
end

% SS1 -> SL1 , SS1 -> PDS & SS1 -> T / SS2 -> SL2 , SS2 -> PDS & SS2 -> T 
% STEP 2 : DETERMINING THE NUMBER OF NUCLEOTIDES ADDED TO SL AND T
if ss(1,end) == "B2c"
    % If length of cs is odd both will execute this statement if they end in 'B2c'
    naSLfSS = [bcount - 1,fcount,(bcount + fcount - 1)]*patchLen;
    naTfSS  = [bcount,fcount,(bcount + fcount - 1)]*patchLen;
else
    % If length of cs is odd both will execute this statement if they end in 'F2c'
    naSLfSS = [fcount,bcount - 1,(bcount + fcount - 1)]*patchLen;
    naTfSS  = [fcount,bcount,(bcount + fcount - 1)]*patchLen;
end
if  y == 1
    extensionTimeArr(1,2) = (naSLfSS/rate);
else
    extensionTimeArr(1,5) = (naSLfSS/rate);
end
% STEP 3 : DETERMINING THE NUMBER OF PDS THAT CAN BE FORMED
primerBindingSites = find((ss == "B2c") | (ss == "F2c"));
primerBindingSites = primerBindingSites(1:(end - 1));
numOfChildPDS      = length(primerBindingSites);
if y == 1
    childPDSnum(1,1) = numOfChildPDS; 
    naChildPDSfSS1   = zeros(0,0); 
else
    childPDSnum(1,3) = numOfChildPDS; 
    naChildPDSfSS2   = zeros(0,0);
end
% STEP 4 : DETERMINING THE COPIES OF PRODUCTS (PDS,SL AND T)
numOfCopiesSL       = round(probSL*numOfCopiesSS);
n                   = numOfChildPDS + 1;
numCopies           = round((numOfCopiesSS - numOfCopiesSL)/n);
if numOfChildPDS > 0
    numOfCopiesPDSc = numCopies;
end
% STEP 5 : IF PDS PRODUCT FORMED, SENDING DATA INTO dataMat AND DATA CONTAINER
if (numOfChildPDS > 0) && (numOfCopiesPDSc > 0)
    for z = 1:numOfChildPDS
        lenOfArrChildPDS = primerBindingSites(z);
        if rem(lenOfArrChildPDS,2) == 0
            bcountPDS   = lenOfArrChildPDS/2;
            fcountPDS   = lenOfArrChildPDS/2;
            naChildPDS  = [bcountPDS,fcountPDS,(bcountPDS + fcountPDS - 1)]*patchLen;
        else
            bcountPDS   = ceil(lenOfArrChildPDS/2); 
            fcountPDS   = lenOfArrChildPDS - bcountPDS;
            if ss(1,lenOfArrChildPDS) == "B2c"
                naChildPDS = [bcountPDS,fcountPDS,(bcountPDS + fcountPDS - 1)]*patchLen;
            else
                naChildPDS = [fcountPDS,bcountPDS,(bcountPDS + fcountPDS - 1)]*patchLen;
            end
        end
        
        ssPDSfSS  = css((end - primerBindingSites(z) + 1):end);
        
        if y == 1
            strPDSfSS1{1,z}      = ssPDSfSS;
            strPDSfSS1{2,z}      = ss; 
            naChildPDSfSS1(1,z)  = naChildPDS;
        else
            strPDSfSS2{1,z}      = ssPDSfSS;
            strPDSfSS2{2,z}      = ss;
            naChildPDSfSS2(1,z)  = naChildPDS;
        end
            
    end
end

% ------------------------------------------------------------------------

% SL1 -> PDS1 / SL2 -> PDS2 
% STEP 1: DETERMINING naPDS
if y == 1
    extensionTimeArr(1,3) = naPDS/rate;
else
    extensionTimeArr(1,6) = naPDS/rate;
end
% STEP 2: DETERMINING NUMBER OF COPIES OF PDS
numOfCopiesPDS = numOfCopiesSL;

% ------------------------------------------------------------------------

% PDS -> SS2, PDS -> PDS, PDS -> SL & PDS -> T / PDS -> SS1, PDS -> PDS,PDS -> SL & PDS -> T
% STEP 1: GENERATING PDS SEQUENCES
ssPDS  = css;
cssPDS = ss;
esPDS  = [cssPDS,ssPDS(2:end)];
lsPDS  = ssPDS(2:end);
revcss = fliplr(cssPDS);
% STEP 2: DETERMINING THE LENGTH OF PDS
lenOfArresPDS = length(esPDS);
% NOTE : lenOfArresPDS will be always odd in this loop 
bcountes  = ceil(lenOfArresPDS/2);
fcountes  = lenOfArresPDS - bcountes;
if esPDS(1,end) == "B2c"
    % If length of ss is odd both the es will execute this statement if they end in 'B2c'
    if y == 1
        naTfPDS1   = [bcountes,fcountes,(bcountes + fcountes - 1)]*patchLen;
        naSLfPDS1  = [bcountes - 1,fcountes,(bcountes + fcountes - 1)]*patchLen;
        lengthPDS1 = naTfPDS1; 
    else
        naTfPDS2   = [bcountes,fcountes,(bcountes + fcountes - 1)]*patchLen;
        naSLfPDS2  = [bcountes - 1,fcountes,(bcountes + fcountes - 1)]*patchLen; 
        lengthPDS2 = naTfPDS2; 
    end
else
    % If length of ss is odd both the es will execute this statement if they end in 'F2c'
    % Then naTfPDS1 = naTfPDS2 and naSLfPDS1 = naSLfPDS2
    if y == 1
        naTfPDS1   = [fcountes,bcountes,(bcountes + fcountes - 1)]*patchLen;
        naSLfPDS1  = [fcountes,bcountes - 1,(bcountes + fcountes - 1)]*patchLen;
        lengthPDS1 = naTfPDS1; 
    else
        naTfPDS2   = [fcountes,bcountes,(bcountes + fcountes - 1)]*patchLen;
        naSLfPDS2  = [fcountes,bcountes - 1,(bcountes + fcountes - 1)]*patchLen;
        lengthPDS2 = naTfPDS2;
    end
end
% STEP 3: DETERMINING THE extensionTimes
if  y == 1
    extensionTimeArr(1,4) = (naSLfPDS1/rate);
else
    extensionTimeArr(1,1) = (naSLfPDS2/rate);
end
% STEP 4: DETERMINING NUMBER OF CHILD PDS 
primerBindingSites        = find((lsPDS == "B2c") | (lsPDS == "F2c"));
primerBindingSites        = primerBindingSites(1:(end - 1));
numOfChildPDS             = length(primerBindingSites);
if y == 1
    childPDSnum(1,2)      = numOfChildPDS;
    naChildPDSfPDS1       = zeros(0,0);
else
    childPDSnum(1,4)      = numOfChildPDS;
    naChildPDSfPDS2       = zeros(0,0);
end
% STEP 5: DETERMINING NUMBER OF COPIES OF PRODUCTS
numOfCopiesSS             = numOfCopiesPDS;
numOfCopiesSL             = round(probSL*numOfCopiesPDS);
n                         = numOfChildPDS + 1;
numCopies                 = round((numOfCopiesPDS - numOfCopiesSL)/n);
if numOfChildPDS > 0
    numOfCopiesPDSc       = numCopies;
end
% STEP 6: SENDING SL AND T INFO INTO dataMat AND DATA CONTAINER
% Sending SL data into data containers
if y == 1
      dsSLfPDS1  = [esPDS;fliplr(ssPDS),fliplr(cssPDS(2:(end - 1))),""];
else
      dsSLfPDS2  = [esPDS;fliplr(ssPDS),fliplr(cssPDS(2:(end - 1))),""];
end
% STEP 7: SENDING CHILD PDS INFO INTO dataMat AND INTO DATA CONTAINER
if (numOfChildPDS > 0) && (numOfCopiesPDSc > 0)
    for z = 1:numOfChildPDS
        lenOfArrChildPDS = lenOfArrSS + primerBindingSites(z);
        if rem(lenOfArrChildPDS,2) == 0
            bcountPDS  = lenOfArrChildPDS/2;
            fcountPDS  = lenOfArrChildPDS/2;
            naChildPDS = [bcountPDS,fcountPDS,(bcountPDS + fcountPDS - 1)]*patchLen; 
        else
            bcountPDS = ceil(lenOfArrChildPDS/2);
            fcountPDS = lenOfArrChildPDS - bcountPDS;
            if esPDS(1,lenOfArrChildPDS) == "B2c"
                naChildPDS = [bcountPDS,fcountPDS,(bcountPDS + fcountPDS - 1)]*patchLen;
            else
                naChildPDS = [fcountPDS,bcountPDS,(bcountPDS + fcountPDS - 1)]*patchLen;
            end 
        end
        
        ssPDSfPDS  = [fliplr(revcss(2:(primerBindingSites(z) + 1))),ssPDS];
        if y == 1
            strPDSfPDS1{1,z}       = ssPDSfPDS;
            strPDSfPDS1{2,z}       = esPDS; 
            naChildPDSfPDS1(1,z)   = naChildPDS;
        else
            strPDSfPDS2{1,z}       = ssPDSfPDS;
            strPDSfPDS2{2,z}       = esPDS;
            naChildPDSfPDS2(1,z)   = naChildPDS;
        end
           
    end
    
end

% --- END OF FOR LOOP 'y' ---
end

% ------------------------------------------------------------------------

% Shorten ssTimeArr to those within simTime

loc = find(ssTimeArr > simTime, 1);
if isempty(loc) == 0
    ssTimeArr   = ssTimeArr(1:(loc - 1));
    ssCopiesArr = ssCopiesArr(1:(loc - 1));
end

dim = size(ssTimeArr);
r   = dim(1,1);

% % Calculating the cycle time
% cycleTime = 0;
% for i = 1:6
%     cycleTime = cycleTime + extensionTimeArr(1,i);
% end
% 
% reqData(1,1) = cycleTime;
% reqData(1,2) = timeToSL1;
% reqData(1,3) = timeToSL2;

% This variable moves from one row to the other
d         = 1;
% These variables are used when dealing with the other SS formation times
shiftCol1 = 0;
shiftCol2 = 0;
shiftCol3 = 0;
shiftCol4 = 0;
% shiftCol is used when you need to move to the next column for the next SS formationTime
shiftCol  = 0;

firstSS = 1;

for w = 1:r
    
while true
    
    if firstSS == 1
        
        numOfCopiesSS         = ssCopiesArr(w,1);
        formationTimeArr(1,1) = ssTimeArr(w,1);
        
        % If the SS coming from slHwy or childPDS has a formation time
        % greater than simTime then break out from loop. If its greater
        % than its parent PDS should come under last amplicon which is
        % taken care of in slHwy and childPDS functions.
        if formationTimeArr(1,1) > simTime
            break;
        end
        
        formationTimeArr(1,2) = formationTimeArr(1,1) + extensionTimeArr(1,2);
        formationTimeArr(1,3) = formationTimeArr(1,2) + extensionTimeArr(1,3);
        formationTimeArr(1,4) = formationTimeArr(1,3) + extensionTimeArr(1,4);
        formationTimeArr(1,5) = formationTimeArr(1,4) + extensionTimeArr(1,5);
        formationTimeArr(1,6) = formationTimeArr(1,5) + extensionTimeArr(1,6);
        
    else
        ssformationTime = round(formationTimeArr(1,1));
        index           = ssformationTime/timeInt + 1;
        if index <= simTime
            dataMat(index,2) = dataMat(index,2) + numOfCopiesSS;
        end
      
        if formationTimeArr(1,1) > simTime
            lastAmpliconsLen(1,la) = 2;
            lastAmpliconsLen(2,la) = lengthPDS2;
            lastAmpliconsLen(3,la) = numOfCopiesPDS;
            la = la + 1;
            break;
        end
    end
    
    % Start executing from else statement in next loop
    firstSS = 0;

    % The number of copies of SS are equal to the number of copies of PDS from where it was formed
    % -------------------------------------------
    numOfCopiesSL  = round(probSL*numOfCopiesSS);
    n              = childPDSnum(1,1) + 1;                        
    numCopies      = round((numOfCopiesSS - numOfCopiesSL)/n);
    if childPDSnum(1,1) > 0
        numOfCopiesPDSc = numCopies;
    end
    numOfCopiesT   = numCopies;
    % -------------------------------------------

    % ADDING CHILD PDS INFO (SS->PDS) INTO dataMat AND INTO DATA CONTAINER
    if (childPDSnum(1,1) > 0) && (numOfCopiesPDSc > 0)
        for i = 1:childPDSnum(1,1)
            childPDSformationTime = round(formationTimeArr(1,1) + (naChildPDSfSS1(1,i)/rate));
            index = childPDSformationTime/timeInt + 1;
            if index <= simTime
                dataMat(index,3) = dataMat(index,3) + numOfCopiesPDSc;
                dataMat(index,6) = dataMat(index,6) + (naChildPDSfSS1(1,i))*numOfCopiesPDSc;
            end
            dataPDSfSS1(d,i + shiftCol1) = formationTimeArr(1,1) + (naChildPDSfSS1(1,i)/rate);
        end
        dataPDSfSS1(d,childPDSnum(1,1) + 1 + shiftCol1) = numOfCopiesPDSc;
    end
    
    % SENDING SL INFO (SS->SL) TO dataMat
    slformationTime = round(formationTimeArr(1,2));
    index = slformationTime/timeInt + 1;
    if index <= simTime
        dataMat(index,4) = dataMat(index,4) + numOfCopiesSL;
        dataMat(index,6) = dataMat(index,6) + extensionTimeArr(1,2)*120*numOfCopiesSL;
    end
    if formationTimeArr(1,2) > simTime
        lastAmpliconsLen(1,la) = 1;
        lastAmpliconsLen(2,la) = lengthSS;
        lastAmpliconsLen(3,la) = numOfCopiesSS;
        la = la + 1;
        break; 
    end

    %SENDING T INFO (SS->T) TO dataMatBlank
    tformationTime = round(formationTimeArr(1,1) + naTfSS/rate);
    index = tformationTime/timeInt + 1;
    if index <= simTime
        dataMat(index,5) = dataMat(index,5) + numOfCopiesT;
        dataMat(index,6) = dataMat(index,6) + naTfSS*numOfCopiesT;
    end
    if ((formationTimeArr(1,1) + naTfSS/rate) <= simTime) && (numOfCopiesT > 0)
        lastAmpliconsLen(1,la) = 4;
        lastAmpliconsLen(2,la) = naTfSS;
        lastAmpliconsLen(3,la) = numOfCopiesT;
        la = la + 1;
    end

    % -------------------------------------------
    numOfCopiesPDS = numOfCopiesSL;
    % -------------------------------------------

    %SENDING PDS INFO (SL->PDS) TO dataMatBlank
    pdsformationTime = round(formationTimeArr(1,3));
    index = pdsformationTime/timeInt + 1;
    if index <= simTime
        dataMat(index,3) = dataMat(index,3) + numOfCopiesPDS;
        dataMat(index,6) = dataMat(index,6) + extensionTimeArr(1,3)*120*numOfCopiesPDS; 
    end
    if formationTimeArr(1,3) > simTime
        lastAmpliconsLen(1,la) = 3;
        lastAmpliconsLen(2,la) = lengthSL;
        lastAmpliconsLen(3,la) = numOfCopiesSL; 
        la = la + 1;
        break;  
    end

    % -------------------------------------------
    numOfCopiesSS = numOfCopiesPDS;
    numOfCopiesSL = round(probSL*numOfCopiesPDS);
    n             = childPDSnum(1,2) + 1;
    numCopies     = round((numOfCopiesPDS - numOfCopiesSL)/n);
    if childPDSnum(1,2) > 0
        numOfCopiesPDSc = numCopies;
    end
    numOfCopiesT = numCopies; 
    % -------------------------------------------

    % ADDING CHILD PDS INFO (PDS->PDS) INTO dataMat AND INTO DATA CONTAINER
    if (childPDSnum(1,2) > 0) && (numOfCopiesPDSc > 0)
        for i = 1:childPDSnum(1,2)
            childPDSformationTime = round(formationTimeArr(1,3) + (naChildPDSfPDS1(1,i)/rate));
            index = childPDSformationTime/timeInt + 1;
            if index <= simTime
                dataMat(index,3) = dataMat(index,3) + numOfCopiesPDSc;
                dataMat(index,6) = dataMat(index,6) + (naChildPDSfPDS1(1,i))*numOfCopiesPDSc;
            end
            dataPDSfPDS1(d,i + shiftCol2) = formationTimeArr(1,3) + (naChildPDSfPDS1(1,i)/rate);
        end
        dataPDSfPDS1(d,childPDSnum(1,2) + 1 + shiftCol2) = numOfCopiesPDSc;
    end

    % ADDING CHILD SL INFO (PDS->SL) INTO SLfPDS1 AND dataMat
    slformationTime = round(formationTimeArr(1,3) + (naSLfPDS1/rate));
    index = slformationTime/timeInt + 1;
    if index <= simTime
        dataMat(index,4) = dataMat(index,4) + numOfCopiesSL;
        dataMat(index,6) = dataMat(index,6) + (naSLfPDS1)*numOfCopiesSL;
    end
    dataSLfPDS1(d,1 + shiftCol) = formationTimeArr(1,3) + (naSLfPDS1/rate);
    dataSLfPDS1(d,2 + shiftCol) = numOfCopiesSL;

    % SENDING SS INFO (PDS->SS) INTO dataMatBlank
    ssformationTime = round(formationTimeArr(1,4));
    index = ssformationTime/timeInt + 1;
    if index <= simTime
        dataMat(index,2) = dataMat(index,2) + numOfCopiesSS;
    end
    if formationTimeArr(1,4) > simTime
        lastAmpliconsLen(1,la) = 2;
        lastAmpliconsLen(2,la) = lengthPDS1;
        lastAmpliconsLen(3,la) = numOfCopiesPDS;
        la = la + 1;
        break;
    end

    % SENDING T INFO (PDS->T) INTO dataMatBlank
    tformationTime = round(formationTimeArr(1,3) + naTfPDS1/rate);
    index = tformationTime/timeInt + 1;
    if index <= simTime
        dataMat(index,5) = dataMat(index,5) + numOfCopiesT;
        dataMat(index,6) = dataMat(index,6) + naTfPDS1*numOfCopiesT;
    end
    if ( (formationTimeArr(1,3) + naTfPDS1/rate) <= simTime) && (numOfCopiesT > 0)
        lastAmpliconsLen(1,la) = 4;
        lastAmpliconsLen(2,la) = naTfPDS1;
        lastAmpliconsLen(3,la) = numOfCopiesT;
        la = la + 1;
    end

    % ------------------------------------------
    numOfCopiesSL  = round(probSL*numOfCopiesSS);
    n              = childPDSnum(1,3) + 1;
    numCopies      = round((numOfCopiesSS - numOfCopiesSL)/n);
    if childPDSnum(1,3) > 0
        numOfCopiesPDSc = numCopies;
    end
    numOfCopiesT   = numCopies;
    % ------------------------------------------

    % ADDING CHILD PDS INFO (SS->PDS) INTO PDSfSS2 AND INTO dataMat
    if (childPDSnum(1,3) > 0) && (numOfCopiesPDSc > 0)
        for i = 1:childPDSnum(1,3)
            childPDSformationTime = round(formationTimeArr(1,4) + (naChildPDSfSS2(1,i)/rate));
            index = childPDSformationTime/timeInt + 1;
            if index <= simTime
                dataMat(index,3) = dataMat(index,3) + numOfCopiesPDSc;
                dataMat(index,6) = dataMat(index,6) + (naChildPDSfSS2(1,i))*numOfCopiesPDSc;
            end
            dataPDSfSS2(d,i + shiftCol3) = formationTimeArr(1,4) + (naChildPDSfSS2(1,i)/rate);
        end
        dataPDSfSS2(d,childPDSnum(1,3) + 1 + shiftCol3) = numOfCopiesPDSc; 
    end

    %SENDING SL INFO (SS->SL) 
    slformationTime = round(formationTimeArr(1,5));
    index = slformationTime/timeInt + 1;
    if index <= simTime
        dataMat(index,4) = dataMat(index,4) + numOfCopiesSL;
        dataMat(index,6) = dataMat(index,6) + extensionTimeArr(1,5)*120*numOfCopiesSL;
    end
    if formationTimeArr(1,5) > simTime
        lastAmpliconsLen(1,la) = 1;
        lastAmpliconsLen(2,la) = lengthSS;
        lastAmpliconsLen(3,la) = numOfCopiesSS;
        la = la + 1;
        break;
    end

    %SENDING T INFO (SS->T) (TERMINATION PRODUCT)
    tformationTime = round(formationTimeArr(1,4) + naTfSS/rate);
    index = tformationTime/timeInt + 1;
    if index <= simTime
        dataMat(index,5) = dataMat(index,5) + numOfCopiesT;
        dataMat(index,6) = dataMat(index,6) + naTfSS*numOfCopiesT;
    end
    if ((formationTimeArr(1,4) + naTfSS/rate) <= simTime) && (numOfCopiesT > 0)
        lastAmpliconsLen(1,la) = 4;
        lastAmpliconsLen(2,la) = naTfSS;
        lastAmpliconsLen(3,la) = numOfCopiesT;
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
        lastAmpliconsLen(1,la) = 3;
        lastAmpliconsLen(2,la) = lengthSL;
        lastAmpliconsLen(3,la) = numOfCopiesSL;
        la = la + 1;
        break;
    end

    % -------------------------------------------
    numOfCopiesSS = numOfCopiesPDS;
    numOfCopiesSL = round(probSL*numOfCopiesPDS);
    n             = childPDSnum(1,4) + 1;
    numCopies     = round((numOfCopiesPDS - numOfCopiesSL)/n);
    if childPDSnum(1,4) > 0
        numOfCopiesPDSc = numCopies;
    end
    numOfCopiesT = numCopies; 
    % -------------------------------------------

    % ADDING CHILD PDS INFO (PDS->PDS) INTO PDSfPDS2
    if (childPDSnum(1,4) > 0) && (numOfCopiesPDSc > 0)
        for i = 1:childPDSnum(1,4)
            childPDSformationTime = round(formationTimeArr(1,6) + (naChildPDSfPDS2(1,i)/rate));
            index = childPDSformationTime/timeInt + 1;
            if index <= simTime  
                dataMat(index,3) = dataMat(index,3) + numOfCopiesPDSc;
                dataMat(index,6) = dataMat(index,6) + (naChildPDSfPDS2(1,i))*numOfCopiesPDSc;
            end
            dataPDSfPDS2(d,i + shiftCol4) = formationTimeArr(1,6) + (naChildPDSfPDS2(1,i)/rate);
        end
        dataPDSfPDS2(d,childPDSnum(1,4) + 1 + shiftCol4) = numOfCopiesPDSc; 
    end

    % ADDING CHILD SL INFO (PDS->SL) INTO SLfPDS2 AND dataMat
    slformationTime      = round(formationTimeArr(1,6) + (naSLfPDS2/rate));
    index                = slformationTime/timeInt + 1;
    if index <= simTime
        dataMat(index,4) = dataMat(index,4) + numOfCopiesSL;
        dataMat(index,6) = dataMat(index,6) + (naSLfPDS2)*numOfCopiesSL;
    end
    dataSLfPDS2(d,1 + shiftCol) = formationTimeArr(1,6) + (naSLfPDS2/rate);
    dataSLfPDS2(d,2 + shiftCol) = numOfCopiesSL;

    %SENDING T INFO (PDS->T) (TERMINATION PRODUCT)
    tformationTime = round(formationTimeArr(1,6) + (naTfPDS2/rate));
    index = tformationTime/timeInt + 1;
    if index <= simTime
        dataMat(index,5) = dataMat(index,5) + numOfCopiesT;
        dataMat(index,6) = dataMat(index,6) + naTfPDS2*numOfCopiesT;
    end
    if ((formationTimeArr(1,6) + (naTfPDS2/rate)) <= simTime) && (numOfCopiesT > 0)
        lastAmpliconsLen(1,la) = 4;
        lastAmpliconsLen(2,la) = naTfPDS2;
        lastAmpliconsLen(3,la) = numOfCopiesT;
        la = la + 1;
    end

    % Updating the formationTimeArr for the next cycle corresponding 
    formationTimeArr(1,1) = formationTimeArr(1,6) + extensionTimeArr(1,1);
    formationTimeArr(1,2) = formationTimeArr(1,1) + extensionTimeArr(1,2);
    formationTimeArr(1,3) = formationTimeArr(1,2) + extensionTimeArr(1,3);
    formationTimeArr(1,4) = formationTimeArr(1,3) + extensionTimeArr(1,4);
    formationTimeArr(1,5) = formationTimeArr(1,4) + extensionTimeArr(1,5);
    formationTimeArr(1,6) = formationTimeArr(1,5) + extensionTimeArr(1,6);

    % Updating the row variable
    d = d + 1;

% --- END OF WHILE LOOP ---
end

% The data will now be entered into the new columns
shiftCol1 = shiftCol1 + childPDSnum(1,1) + 1;
shiftCol2 = shiftCol2 + childPDSnum(1,2) + 1;
shiftCol3 = shiftCol3 + childPDSnum(1,3) + 1;
shiftCol4 = shiftCol4 + childPDSnum(1,4) + 1;
shiftCol  = shiftCol  + 2; 
d = 1;

firstSS = 1;

end

% This statement will be executed if in the loop you have not covered 6 amplicons OR a break from for loop 'w'

if breakTheLoop == 1
        dataFun{1}  = strPDSfSS1;
        dataFun{2}  = dataPDSfSS1;
        dataFun{3}  = strPDSfPDS1;
        dataFun{4}  = dataPDSfPDS1;
        dataFun{5}  = strPDSfSS2;
        dataFun{6}  = dataPDSfSS2;
        dataFun{7}  = strPDSfPDS2;
        dataFun{8}  = dataPDSfPDS2;
        dataFun{9}  = dsSLfPDS1;
        dataFun{10} = dsSLfPDS2;
        dataFun{11} = dataSLfPDS1;
        dataFun{12} = dataSLfPDS2;
        dataFun{13} = childPDSnum;
        dataFun{14} = lastAmpliconsLen;
        dataFun{15} = dataMat;
        return;
else
    if x == 1
        dataFun{1}  = strPDSfSS1;
        dataFun{2}  = dataPDSfSS1;
        dataFun{3}  = strPDSfPDS1;
        dataFun{4}  = dataPDSfPDS1;
        dataFun{5}  = strPDSfSS2;
        dataFun{6}  = dataPDSfSS2;
        dataFun{7}  = strPDSfPDS2;
        dataFun{8}  = dataPDSfPDS2;
        dataFun{9}  = dsSLfPDS1;
        dataFun{10} = dsSLfPDS2;
        dataFun{11} = dataSLfPDS1;
        dataFun{12} = dataSLfPDS2;
        dataFun{13} = childPDSnum;
        dataFun{14} = lastAmpliconsLen;
        dataFun{15} = dataMat;
        % Go to the next iteration in for loop x
        continue;
    else
        dataFun{16}  = strPDSfSS1;
        dataFun{17}  = dataPDSfSS1;
        dataFun{18}  = strPDSfPDS1;
        dataFun{19}  = dataPDSfPDS1;
        dataFun{20}  = strPDSfSS2;
        dataFun{21}  = dataPDSfSS2;
        dataFun{22}  = strPDSfPDS2;
        dataFun{23}  = dataPDSfPDS2;
        dataFun{24}  = dsSLfPDS1;
        dataFun{25}  = dsSLfPDS2;
        dataFun{26}  = dataSLfPDS1;
        dataFun{27}  = dataSLfPDS2;
        dataFun{28}  = childPDSnum;
        dataFun{29}  = lastAmpliconsLen;
        dataFun{30}  = dataMat;
        % Break out of for loop x
        return
    end       
end     

% --- END OF FOR LOOP WHICH CONSIDERS BOTH THE STRANDS ---
end

% --- END OF FUNCTION ---
end