function MS = msbinoc(MSleft,MSright,binT)

% binT: limit of the time difference between two MS onsets (left/right)
%%
param = fieldnames(MSleft);

if ~isempty(MSleft.onset) && ~isempty(MSright.onset)
    cMat = abs(repmat(MSleft.onset,1,length(MSright.onset)) - repmat(MSright.onset',length(MSleft.onset),1));
    if size(cMat,1)>1 && size(cMat,2)>1
    [lMat,leftindex] = min(cMat);
    [rMat,rightindex] = min(cMat'); %#ok<UDIM>
    elseif size(cMat,1)==1 && size(cMat,2)>1
        lMat = cMat';
        leftindex = ones(size(lMat));
        [rMat,rightindex] = min(cMat'); %#ok<UDIM>
    elseif size(cMat,1)>1 && size(cMat,2)==1
        [lMat,leftindex] = min(cMat);
        rMat = cMat';
        rightindex = ones(size(rMat));
    elseif size(cMat,1)==1 && size(cMat,2)==1
        lMat = cMat;
        rMat = cMat;
        leftindex = 1;
        rightindex = 1;
    end
    lx = find(lMat<binT);
    rx = find(rMat<binT);
    L = leftindex(lx);
    R = rightindex(rx);
    LMat = lMat(lx);
    RMat = rMat(rx);

    remL=[];
    remR=[];
    if length(L)~=length(R)
        x = find(diff(L)==0);
        if isempty(x)==0
            remL = zeros(1,length(x));
            for k=1:length(x)
                [~,b] = max([LMat(x(k)) LMat(x(k)+1)]);
                remL(k) = x(k)+b-1;
            end
            L(remL)=[]; %#ok<*NASGU>
        end
        x = find(diff(R)==0);
        if isempty(x)==0
            remR = zeros(1,length(x));
            for k=1:length(x)
                [~,b] = max([RMat(x(k)) RMat(x(k)+1)]);
                remR(k) = x(k)+b-1;
            end
            R(remR)=[];
        end
    end
    % 20151217下記追加。researchnotebookno.3,p.28-31に詳細記載。
    % 削除対象となったグループに対して、反対側でそのグループのいずれかの要素を最近傍としていた場合には、それを削除しないといけない。
    cRx = common(R,lx(remL));
    cLx = common(L,rx(remR));
    for k=1:length(cRx)
        R(R==cRx(k))=[];
    end
    for k=1:length(cLx)
        L(L==cLx(k))=[];
    end
    while length(L)~=length(R)
        tmp1 = [];
        for k=1:length(cRx)
            tmp1 = [tmp1 find(rightindex==cRx(k))]; %#ok<AGROW>
        end
        tmp2 = [];
        for k=1:length(cLx)
            tmp2 = [tmp2 find(leftindex==cLx(k))]; %#ok<AGROW>
        end
        for k=1:length(tmp1)
            L(L==tmp1(k))=[];
        end
        for k=1:length(tmp2)
            R(R==tmp2(k))=[];
        end
        cLx = tmp1;
        cRx = tmp2;
    end

    for k=1:length(param)
        eval(['MS.binoc.',param{k},' = (MSleft.',param{k},'(L) + MSright.',param{k},'(R))/2;'])
        eval(['MS.left.',param{k},' = MSleft.',param{k},'(L);'])
        eval(['MS.right.',param{k},' = MSright.',param{k},'(R);'])
    end

else
    for k=1:length(param)
        eval(['MS.binoc.',param{k},' = [];'])
        eval(['MS.left.',param{k},' = [];'])
        eval(['MS.right.',param{k},' = [];'])
    end
end
