function [dataFun] = slHwy(slPackage1,slPackage2,dataMat,simTime,timeInt,rate,probSL,patchLen)

% Outline of this function is 
% SL1 -> PDS1 -> SL2 -> PDS2 ...
% 1. This highway continues producing larger and larger SL and PDS.
% 2. This is an an important aspect as SL1 is not produced again by any mechanism in this pathway as it is a highway. Hence we can consider the entire
%    sltimeArr at once and not take individual elements as in ssloop function. In ssloop SS1 is being generated back after a fixed period of time
% 3. OUTPUTS OF THIS FUNCTION ARE SS AND PDS  

% Outer for loop takes care of the sister amplicons formed from pathway 1 and pathway 2
% EX : [B2  F2  B2c ] and [F2  B2  F2c]
%      [B2c F2c     ]     [F2c B2c    ]

% CASE 1 : NO AMPLICONS ARE INPUT
if (isempty(slPackage1) && isempty(slPackage2))
    dataFun = cell(0,0);
    return
% CASE 2 : AMPLICON FROM PATHWAY 1 IS INPUT
elseif (isempty(slPackage1))
    ds            = slPackage2{1};
    slTimeArr     = slPackage2{2}; 
    slCopiesArr   = slPackage2{3};
    breakTheLoop  = 1;
% CASE 3 : AMPLICON FROM PATHWAY 2 IS INPUT
elseif (isempty(slPackage2))
    ds            = slPackage1{1};
    slTimeArr     = slPackage1{2}; 
    slCopiesArr   = slPackage1{3};
    breakTheLoop  = 1;
% CASE 4 : BOTH AMPLICONS ARE INPUT
else
    breakTheLoop  = 0; 
end

% breakTheLoop variable checks how many amplicon inputs are provided to fun 

% Data variables which store data of both sister strands
% Number coding used to identify the type of amplicon 
% 1 - SS
% 2 - PDS
% 3 - SL
% 4 - T
lastAmpliconsLen = zeros(0,0);
la = 1;

for x = 1:2

if breakTheLoop == 0
    if  x == 1
        ds            = slPackage1{1};
        slTimeArr     = slPackage1{2};
        slCopiesArr   = slPackage1{3};
    else
        ds            = slPackage2{1};
        slTimeArr     = slPackage2{2}; 
        slCopiesArr   = slPackage2{3};
    end
end

% Creating the cells and arrays used in this function
%     -- SSfPDS  --
strSSfPDS   = cell(0,0);
dataSSfPDS  = zeros(0,0);
%     -- PDSfPDS --
strPDSfPDS  = cell(0,0);
dataPDSfPDS = zeros(0,0);
% The below array will contain the number of child PDS that can be produced from a parent PDS
childPDSnum = zeros(0,0);

% Variables used in these data containers
shiftCol1 = 0;
shiftCol2 = 0;
% shiftCol3 used in dataSSfPDS
shiftCol3 = 0;
loadSS    = 1;
cpn       = 1;

% This variable is used as a switch telling the program where to get its SL timeStamps from
isFirst = 1;

% First remove any '0' values from the slTimeArray and slCopyArray
loc = find(slTimeArr == 0,1);
if isempty(loc) == 0
    slCopiesArr = slCopiesArr(1:(loc - 1));
    slTimeArr   = slTimeArr(1:(loc - 1));
end

r = length(slTimeArr);

% This while loop keeps running till a larger PDS cannot be formed
while true
    
% SL -> PDS
% STEP 1 : CHECK IF FIRST ENTRY IS LESS THAN SIMTIME
if slTimeArr(1,1) > simTime
    break;
end
% STEP 2 : OBTAINING THE SL FORMATION TIME
% If its first cycle then you get it from function input otherwise from pdsformationTime
if isFirst == 0
    % Determining the next SL sequence (VERIFIED)
    ds = [esPDS;fliplr(ssPDS),fliplr(cssPDS(2:(end - 1))),""];
    slTimeArr = pdsTimeArr + (naSL/rate);
    for y = 1:r
        slTime = round(slTimeArr(y,1));
        index  = slTime/timeInt + 1;
        if index <= simTime
            dataMat(index,4) = dataMat(index,4) + slCopiesArr(y,1);
            dataMat(index,6) = dataMat(index,6) + naSL*slCopiesArr(y,1);
        else
            break;
        end    
    end
end
% Changing the counter value so that from next iteration SL can now add data into dataMat
isFirst = 0;
% STEP 3 : DETERMINING THE LENGTH OF SL AND THE NUCLEOTIDES ADDED TO PDS
% The ds structures always have an odd length
lenOfArrSL = length(ds);
bcount     = ceil(lenOfArrSL/2);
fcount     = lenOfArrSL - bcount;
if ds(1,end) == "B2c"
    lengthSL  = [bcount,fcount,(bcount + fcount - 1)]*patchLen;
    naPDS     = lengthSL;
else
    lengthSL  = [fcount,bcount,(bcount + fcount - 1)]*patchLen;
    naPDS     = lengthSL; 
end
% STEP 4 : DETERMINING THE PDS FORMATION TIME
pdsTimeArr       = slTimeArr + (naPDS/rate);
extensionTimePDS = (naPDS/rate);
% STEP 5 : DETERMINING THE NUMBER COPIES OF PDS
pdsCopiesArr     = slCopiesArr;

% PDS -> SS, PDS -> PDS, PDS -> SL and PDS -> T
% STEP 1 : SENDING PDS INFO INTO dataMat 
for y = 1:r
    
    pdsTime = round(pdsTimeArr(y,1));
    index   = pdsTime/timeInt + 1;
    
    if index <= simTime
        dataMat(index,3) = dataMat(index,3) + pdsCopiesArr(y,1);
        dataMat(index,6) = dataMat(index,6) + naPDS*pdsCopiesArr(y,1);
    end
    
    if pdsTimeArr(y,1) > simTime
        timeSLwasFormed = pdsTimeArr(y,1) - extensionTimePDS;
        if timeSLwasFormed <= simTime
            lastAmpliconsLen(1,la) = 3;
            lastAmpliconsLen(2,la) = lengthSL;
            lastAmpliconsLen(3,la) = slCopiesArr(y,1);
            la = la + 1;
        end
    end
    
end
% If first element is greater than simTime break from loop
if pdsTimeArr(1,1) > simTime
    break;
end
% STEP 2 : GENERATE THE PDS SEQUENCE AND FIND ITS PRIMER BINDING SITES
if ds(1,end) == "B2c"
ssPDS  = ["B2",fliplr(ds(2,1:(end-1)))];
lengthOfSS = length(ssPDS);
esPDS  = [ds(1,1:end),fliplr(ds(2,1:(end-1)))];
cssPDS = esPDS(1:lengthOfSS);
revcss = fliplr(cssPDS);
lsPDS  = fliplr(ds(2,1:(end-1)));
else
ssPDS  = ["F2",fliplr(ds(2,1:(end-1)))];
lengthOfSS = length(ssPDS);
esPDS  = [ds(1,1:end),fliplr(ds(2,1:(end-1)))];
cssPDS = esPDS(1:lengthOfSS);
revcss = fliplr(cssPDS);
lsPDS  = fliplr(ds(2,1:(end-1)));
end
% STEP 3 : DETERMINING naSL,naT and lengthPDS
lenOfArresPDS = length(esPDS);
bcountes      = ceil(lenOfArresPDS/2);
fcountes      = lenOfArresPDS - bcountes;
if esPDS(1,end) == "B2c"
    naSL = [(bcountes - 1),fcountes,(bcountes + fcountes - 1)]*patchLen;
    naT  = [bcountes,fcountes,(bcountes + fcountes - 1)]*patchLen; 
    lengthPDS = naT;
else
    naSL = [fcountes,(bcountes - 1),(bcountes + fcountes - 1)]*patchLen;
    naT  = [fcountes,bcountes,(bcountes + fcountes - 1)]*patchLen; 
    lengthPDS = naT;
end
% STEP 4 : DETERMINING ssFormationTime
% NOTE : We take ssformationTime as the strand displacement time to form SL
ssTimeArr        = pdsTimeArr + (naSL/rate);
tTimeArr         = pdsTimeArr + (naT/rate);
displacementTime = (naSL/rate);
% STEP 5 : Searching for primer binding sites
primerBindingSites = find((lsPDS == "B2c") | (lsPDS == "F2c"));
primerBindingSites = primerBindingSites(1:(end - 1));
numOfChildPDS      = length(primerBindingSites);
% STEP 6 : Determining the number of copies of SS,PDS,SL and T that could be formed from parent PDS
ssCopiesArr = pdsCopiesArr;
slCopiesArr = round(probSL*pdsCopiesArr);
n           = numOfChildPDS + 1;
numCopies   = round((pdsCopiesArr - slCopiesArr)/n);
if numOfChildPDS > 0
    pdsChildCopiesArr = numCopies;
end
tCopiesArr  = numCopies;
% Can enter the loop only if num of copies of child PDS > 0
if (numOfChildPDS > 0) && (pdsChildCopiesArr(1,1) > 0)
    for y = 1:numOfChildPDS
        lenOfArrChildPDS = lenOfArrSL + primerBindingSites(y);
        if rem(lenOfArrChildPDS,2) == 0
            bcountPDS    = lenOfArrChildPDS/2;
            fcountPDS    = lenOfArrChildPDS/2;
            naChildPDS   = [bcountPDS,fcountPDS,(bcountPDS + fcountPDS - 1)]*patchLen;
        else
            bcountPDS    = ceil(lenOfArrChildPDS/2);
            fcountPDS    = lenOfArrChildPDS - bcountPDS;
            if esPDS(1,lenOfArrChildPDS) == "B2c"
                naChildPDS = [bcountPDS,fcountPDS,(bcountPDS + fcountPDS - 1)]*patchLen;
            else
                naChildPDS = [fcountPDS,bcountPDS,(bcountPDS + fcountPDS - 1)]*patchLen;
            end 
        end

        childPDSTimeArr = pdsTimeArr + (naChildPDS/rate);

        % SENDING CHILD PDS INFO INTO dataMat AND INTO DATA CONTAINER         
        % Sending into dataMat
        for z = 1:r
            childpdsTime = round(childPDSTimeArr(z,1));
            index        = childpdsTime/timeInt + 1;
            if index <= simTime
                dataMat(index,3) = dataMat(index,3) + pdsChildCopiesArr(z,1);
                dataMat(index,6) = dataMat(index,6) + naChildPDS*pdsChildCopiesArr(z,1);
            else
                break;
            end
        end

        % Sending child PDS information into data containers
        ssPDSfPDS                        = [fliplr(revcss(2:(primerBindingSites(y) + 1))),ssPDS];
        strPDSfPDS{1,y + shiftCol2}      = ssPDSfPDS;
        strPDSfPDS{2,y + shiftCol2}      = esPDS;
        dataPDSfPDS(1:r,y + shiftCol1)   = childPDSTimeArr;               
    end
    childPDSnum(1,cpn)                             = numOfChildPDS;
    dataPDSfPDS(1:r,numOfChildPDS + 1 + shiftCol1) = pdsChildCopiesArr;
    shiftCol1 = shiftCol1 + numOfChildPDS + 1;
    shiftCol2 = shiftCol2 + numOfChildPDS;
    cpn = cpn + 1;
end
% STEP 7 : 6.1 SENDING CHILD SS INFO INTO datMat AND INTO DATA CONTAINER  
% Sending into dataMat

for y = 1:r
    
    ssTime = round(ssTimeArr(y,1));
    index  = ssTime/timeInt + 1;
    
    if index <= simTime
        dataMat(index,2) = dataMat(index,2) + ssCopiesArr(y,1);
    end
    
    if ssTimeArr(y,1) > simTime
        timePDSwasFormed = ssTimeArr(y,1) - displacementTime;
        if timePDSwasFormed <= simTime
            lastAmpliconsLen(1,la) = 2;
            lastAmpliconsLen(2,la) = lengthPDS;
            lastAmpliconsLen(3,la) = pdsCopiesArr(y,1);
            la = la + 1;
        end
    end
    
end

% Sending into data container
strSSfPDS{1,loadSS}                 = ssPDS;
strSSfPDS{2,loadSS}                 = cssPDS;
dataSSfPDS(1:r,1 + shiftCol3)       = ssTimeArr;
dataSSfPDS(1:r,2 + shiftCol3)       = ssCopiesArr;
shiftCol3 = shiftCol3 + 2;
loadSS    = loadSS + 1;

% STEP 7 : SENDING T INFO INTO dataMat

for y = 1:r 
    
    tTime = round(tTimeArr(y,1));
    index = tTime/timeInt + 1;
    
    if index <= simTime
        dataMat(index,5) = dataMat(index,5) + tCopiesArr(y,1);
        dataMat(index,6) = dataMat(index,6) + naT*tCopiesArr(y,1);
    end
    
    if (tTimeArr(y,1) <= simTime) && (tCopiesArr(y,1) > 0)
        lastAmpliconsLen(1,la) = 4;
        lastAmpliconsLen(2,la) = naT;
        lastAmpliconsLen(3,la) = tCopiesArr(y,1);
        la = la + 1;
    end
    
end

% ---WHILE LOOP END---
end

% If a break occurs in the while loop it jumps here

if breakTheLoop == 1
    dataFun{1} = strSSfPDS;
    dataFun{2} = dataSSfPDS;
    dataFun{3} = strPDSfPDS;
    dataFun{4} = dataPDSfPDS;
    dataFun{5} = childPDSnum;
    dataFun{6} = lastAmpliconsLen;
    dataFun{7} = dataMat;
    break;
else
% If it breaks out of the while loop and lands here then you have to add the data containers to dataFun and pass control to next iteration in for loop 'x'
    if x == 1
        dataFun{1}  = strSSfPDS;
        dataFun{2}  = dataSSfPDS;
        dataFun{3}  = strPDSfPDS;
        dataFun{4}  = dataPDSfPDS;
        dataFun{5}  = childPDSnum;
        dataFun{6}  = lastAmpliconsLen;
        dataFun{7}  = dataMat;
        % Passing control to for loop 'x' for next iteration
        continue;
    else
        dataFun{8}   = strSSfPDS;
        dataFun{9}   = dataSSfPDS;
        dataFun{10}  = strPDSfPDS;
        dataFun{11}  = dataPDSfPDS;
        dataFun{12}  = childPDSnum;
        dataFun{13}  = lastAmpliconsLen;
        dataFun{14}  = dataMat;
        % Jump out of for loop 'x'
        break;
    end
end

% ---FOR LOOP 'x' END  ---
end

% ---FUNCTION END---
end

