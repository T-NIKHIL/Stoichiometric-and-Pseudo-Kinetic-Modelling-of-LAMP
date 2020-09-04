function [dataFun] = childPDS(pdsPackage1,pdsPackage2,dataMat,simTime,timeInt,rate,probSL,patchLen)

% Outline of this function is 
% PDSn -> SS1 , PDSn -> SL , PDSn -> T
% 1. Here 'n' represents the child PDS number. The child PDS can be formed from SS or PDS 
% 2. Here we assume that the child PDS's smaller than the largest child PDS do not form the remaining larger child PDS
% 3. The above assumption is valid as we take very larger probabilities for formation of SL. Hence almost all of the reactant goes to form SL
% 4. OUTPUTS OF THIS FUNCTION ARE SS AND SL

% Need to first check if any of the packages are empty
% CASE 1 : NO AMPLICONS ARE INPUT
if (isempty(pdsPackage1) && isempty(pdsPackage2))
    dataFun = cell(0,0);
    return
% CASE 2 : AMPLICON FROM PATHWAY 1 IS INPUT
elseif (isempty(pdsPackage1))
    % If pdsPackage1 is empty then take 
    ssPDS        = pdsPackage2{1};
    esPDS        = pdsPackage2{2};
    pdsTimeArr   = pdsPackage2{3};
    pdsCopiesArr = pdsPackage2{4};
    breakTheLoop = 1;
% CASE 3 : AMPLICON FROM PATHWAY 2 IS INPUT
elseif (isempty(pdsPackage2))
    ssPDS        = pdsPackage1{1};
    esPDS        = pdsPackage1{2};
    pdsTimeArr   = pdsPackage1{3};
    pdsCopiesArr = pdsPackage1{4};
    breakTheLoop = 1;
% CASE 4 : BOTH AMPLICONS ARE INPUT
else
    breakTheLoop = 0;
end

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
    if x == 1
        ssPDS        = pdsPackage1{1};
        esPDS        = pdsPackage1{2};
        pdsTimeArr   = pdsPackage1{3};
        pdsCopiesArr = pdsPackage1{4};
    else
        ssPDS        = pdsPackage2{1};
        esPDS        = pdsPackage2{2};
        pdsTimeArr   = pdsPackage2{3};
        pdsCopiesArr = pdsPackage2{4};
    end
end

% Creating the data containers
%     -- SSfPDS  --
strSSfPDS  = cell(0,0);
dataSSfPDS = zeros(0,0);
%     -- SLfPDS  -- 
dataSLfPDS = zeros(0,0);

% STEP 1 : SHORTENING THE ARRAY TO CUT OUT THE PORTION CONTAINING 0 COPIES OF PDS
loc = find(pdsCopiesArr == 0,1);
if isempty(loc) == 0
    pdsCopiesArr = pdsCopiesArr(1:(loc - 1));
    pdsTimeArr   = pdsTimeArr(1:(loc - 1));
end

r = length(pdsTimeArr);

% STEP 3 : GENERATING THE PDS SEQUENCES 
lengthOfSS     = length(ssPDS);
lengthOfES     = length(esPDS);
lengthOfCSS    = lengthOfES - lengthOfSS;
% This is required for genearting the SL structure from PDS
cssPDS         = esPDS(1:lengthOfCSS);
% Required for checking the number of PDS that could be formed
lsPDS          = esPDS((lengthOfSS + 1):end);
% This is for the SS produced from PDS
ss             = ssPDS;
css            = esPDS(1:lengthOfSS);

% STEP 4 : DETERMINING lengthPDS, naSL AND naT
if rem(lengthOfES,2) == 0
    bcountes  = lengthOfES/2;
    fcountes  = lengthOfES/2;
    naT       = [bcountes,fcountes,(bcountes + fcountes - 1)]*patchLen;
    lengthPDS = naT; 
    if esPDS(1,end) == "B2c"
        naSL = [(bcountes - 1),fcountes,(bcountes + fcountes - 1)]*patchLen;
    else
        naSL = [bcountes,(fcountes - 1),(bcountes + fcountes - 1)]*patchLen;
    end
else
    bcountes = ceil(lengthOfES/2);
    fcountes = lengthOfES - bcountes;
    naSL     = [(bcountes - 1),fcountes,(bcountes + fcountes - 1)]*patchLen;
    if esPDS(1,end) == "B2c"
        naT  = [bcountes,fcountes,(bcountes + fcountes - 1)]*patchLen;
        lengthPDS = naT;
    else
        naT  = [fcountes,bcountes,(bcountes + fcountes - 1)]*patchLen;
        lengthPDS = naT;
    end
end
% STEP 5 : DETERMINING ssformationTime, slformatonTime, tformationTime
% HERE WE CONSIDER THE STRAND DISPLACEMENT TIME CORRESPONDING TO SL FOR THE FORMATION OF SS
ssformationTime  = pdsTimeArr + (naSL/rate);
displacementTime = (naSL/rate); 
slformationTime  = pdsTimeArr + (naSL/rate);
tformationTime   = pdsTimeArr + (naT/rate);
% STEP 6 : CREATING THE STRUCTURE OF SL
ds = [esPDS;fliplr(ssPDS),fliplr(cssPDS(2:end)),""];
% We already have the SS structure which will get produced from PDS. It is nothing but ssPDS
% STEP 7 : FINDING THE PRIMER BINDING SITES
primerBindingSites = find((lsPDS == "B2c") | (lsPDS == "F2c"));
primerBindingSites = primerBindingSites(1:(end - 1));
numOfChildPDS      = length(primerBindingSites);
% STEP 8 : DETERMINING THE COPIES OF PRODUCTS FORMED
ssCopiesArr    = pdsCopiesArr;
slCopiesArr    = round(probSL*pdsCopiesArr);
n              = numOfChildPDS + 1;
numCopies      = round((pdsCopiesArr - slCopiesArr)/n);
tCopiesArr     = numCopies;
% STEP 9 : SENDING T INTO dataMat
for y = 1:r
    tTime = round(tformationTime(y,1));
    index = tTime/timeInt + 1;
    if index <= simTime
        dataMat(index,5) = dataMat(index,5) + tCopiesArr(y,1);
        dataMat(index,6) = dataMat(index,6) + naT*tCopiesArr(y,1);
    end
    if (tTime <= simTime) && (tCopiesArr(y,1) > 0)
        lastAmpliconsLen(1,la) = 4;
        lastAmpliconsLen(2,la) = naT;
        lastAmpliconsLen(3,la) = tCopiesArr(y,1);
        la = la + 1;
    end
end

% STEP 10 : SENDING SS AND SL INTO dataMat AND INTO DATA CONTAINERS

% Sending data into dataMat

for y = 1:r
    ssTime = round(ssformationTime(y,1));
    index  = ssTime/timeInt + 1;
    if index <= simTime
        dataMat(index,2) = dataMat(index,2) + ssCopiesArr(y,1);
    end
    if ssformationTime(y,1) > simTime
        timePDSwasFormed = ssformationTime(y,1) - displacementTime;
        if timePDSwasFormed <= simTime
            lastAmpliconsLen(1,la) = 2;
            lastAmpliconsLen(2,la) = lengthPDS;
            lastAmpliconsLen(3,la) = pdsCopiesArr(y,1);
            la = la + 1;
        end  
    end
end

for y = 1:r
    slTime = round(slformationTime(y,1));
    index  = slTime/timeInt + 1;
    if index <= simTime
        dataMat(index,4) = dataMat(index,4) + slCopiesArr(y,1);
        dataMat(index,6) = dataMat(index,6) + naSL*slCopiesArr(y,1);
    else
        break;
    end
end

% Sending data into data containers

strSSfPDS{1,1}           = ss;
strSSfPDS{2,1}           = css;
dataSSfPDS(1:r,1)        = ssformationTime;
dataSSfPDS(1:r,2)        = ssCopiesArr;
dsSLfPDS                 = ds;
dataSLfPDS(1:r,1)        = slformationTime;
dataSLfPDS(1:r,2)        = slCopiesArr;

% STEP 11 : SENDING DATA INTO dataFun CONTAINER
if breakTheLoop == 1
    dataFun{1} = strSSfPDS;
    dataFun{2} = dataSSfPDS;
    dataFun{3} = dsSLfPDS;
    dataFun{4} = dataSLfPDS;
    dataFun{5} = lastAmpliconsLen;
    dataFun{6} = dataMat;
    break;
else
    if x == 1
        dataFun{1}  = strSSfPDS;
        dataFun{2}  = dataSSfPDS;
        dataFun{3}  = dsSLfPDS;
        dataFun{4}  = dataSLfPDS;
        dataFun{5}  = lastAmpliconsLen;
        dataFun{6}  = dataMat;
        % Continue to next iteration
        continue;
    else
        dataFun{7}  = strSSfPDS;
        dataFun{8}  = dataSSfPDS;
        dataFun{9}  = dsSLfPDS;
        dataFun{10} = dataSLfPDS;
        dataFun{11} = lastAmpliconsLen;
        dataFun{12} = dataMat;
        break;
    end
end

% END OF FOR LOOP
end

% END OF FUNCTION
end

