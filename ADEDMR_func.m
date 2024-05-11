function [gbest,gbestval, nfes,RecordT] = ADEDMR_func(fhd,D, max_nfes, Xmin, Xmax, varargin)
% Parameter Adaptive Learning Mechanism DE (palm-DE)
% F = 0.5, CR = 0.5, numst = D;
%%
% random seed
%rand('state',sum(100*clock));
stm = RandStream('swb2712','Seed',sum(100*clock));
RandStream.setGlobalStream(stm);

%argument passing
%     F = 0.5;
%     Cr = 0.5;
ps_ini = 18*D;
%     ps_ini = round(25*log(D)*sqrt(D));
ps_min = 10;
ps = ps_ini;

% if length(Xmin) == 1
%     Rmin = repmat(Xmin,1,D);
%     Rmax = repmat(Xmax,1,D);
% end
% Xmin= repmat(Rmin,ps,1);
% Xmax = repmat(Rmax,ps,1);

fidvec = cell2mat(varargin);
fid = fidvec(1);
runid = fidvec(2);
WriteFlag = true;
if WriteFlag
    targetbest = [100;200;300;400;500;600;700;800;900;1000;1100;1200;1300;1400;1500;1600;1700;1800;1900;
        2000;2100;2200;2300;2400;2500;2600;2700;2800;2900;3000];
    name = ['ADEDMR_fid_',num2str(fid),'_',num2str(D),'D_',num2str(runid),'.dat'];
    fout = fopen(name,'a');
end
tic;
%initialize;
pos = Xmin+ (Xmax-Xmin).*rand(ps,D);
eval = feval(fhd,pos',varargin{:});
nfes = ps;
pbest = pos;
pbestval = eval;
[gbestval,gbestid] = min(pbestval);
[gworstval,~] = max(pbestval);
gbest = pos(gbestid,:);
if WriteFlag
    fprintf(fout,'%d\t%.15f\n',1,gbestval-targetbest(fid));
end

interval_num=10000*D/20;
print_num_record = interval_num;
Archfactor = 1.6;
% ArchfactorB_max = 0.6;
% ArchfactorB_min = 0.2 ;
ArchfactorB = 0.5;
pbest_rate_max = 0.2;
pbest_rate_min = 0.05;
pbest_rate =pbest_rate_max;
memory_size = 5;
memory_MUF = 0.5.*ones(memory_size,1);
memory_MUCr = 0.5.*ones(memory_size,1);
memory_order = 1;
decayA = [];
% T0 = 1.0;
% pmid = zeros(ps,D);
% decayRate = T0/70;
A = [];
B = [];
g = 1;
%     pivot = [0/3*max_nfes,3/3*ps_ini];
Xmax = 100;
Xmin = -100;
counter = zeros(ps,1);
% Icount = zeros(ps,1);
%     Mean1 = mean(pbest);
n = 40;
% dcount = 0;
%     RDI=1;
%     DI_ini = sqrt(sum(pdist2(Mean1,pbest)));
vlim = log(1+200.^D);

%     Seedval_store = [];
probStra = zeros(memory_size,1);


lu=[-100*ones(1,D);100*ones(1,D)];

% for i=1:ps
%     neighbindex=mod([ps+i-1-N_rsize(i):ps+i-1+N_rsize(i)],ps)+1;
%     [nfit,nind]=min(pbestval(neighbindex));
%     Gfit=[Gfit nfit];
%     Gind=[Gind neighbindex(nind)];
% end

ps_WW = 10;
% popLS = Xmin+ (Xmax-Xmin).*rand(ps_WW,D);
% popLSeval = feval(fhd,popLS',varargin{:});
% [~,indexLS] = sort(popLSeval);
% BestPoint = popLS(indexLS(1), :);
% nfes = nfes + ps_WW;
%         %%% Initialize LS population for re-start
% mpbest_rate_min = 0.5;
% mpbest_rate_max = 0.8;
%  mpbest_rate = 0.5;
% pmid2 = [];
[~,indexSel] = sort(pbestval);
while nfes < max_nfes
    %          Distance = pdist2(gbest,pbest);
    Fd=[];
    pps = max(round(pbest_rate*ps),2);
    lenSel = max(ceil(pps*rand(1,ps)),1);
    pbestIndex = zeros(ps,1);
    for idx = 1:ps
        pbestIndex(idx) = indexSel(lenSel(idx));
    end
    pbestB = pos(pbestIndex,:);
    
    
    %         mps = round(((1-mpbest_rate)/2)*ps) + max(round(mpbest_rate*ps),2);
    %         lenSel = max(ceil(mps*rand(1,ps)),1);
    %         pbestIndex = zeros(ps,1);
    %         for idx = 1:ps
    %             pbestIndex(idx) = indexSel(lenSel(idx));
    %         end
    %         pbestM = pos(pbestIndex,:);
    
    
    %     pmid(indexSel(ps+1:length(indexSel))',:) = [];
    %     pmidval(indexSel(ps+1:length(indexSel))') = [];
    %     pmidval_id(indexSel(ps+1:length(indexSel))') = [];
    
    
    popLS = pos(indexSel(1:ps_WW),:);
    BestPoint = pos(indexSel(1), :);
    
    %%
    rndBase = randperm(ps)';
    %       rndBase1 = randperm(ps)';
    psExt = ps + size(A,1);
    rndSeq1 = ceil(rand(ps,1)*psExt);
    %     psExt1 = ps + size(pmid2,1);
    %     rndSeq2 = ceil(rand(ps,1)*psExt1);
    psExtB = ps + size(B,1);
    rndSeq2 = ceil(rand(ps,1)*psExtB);
    for ii = 1:ps
        while rndBase(ii)==ii
            rndBase(ii)=ceil(rand()*ps);
        end
        %                     while rndBase1(ii)==ii || rndBase1(ii)== rndBase(ii)
        %                         rndBase1(ii)=ceil(rand()*ps);
        %                     end
        while  rndSeq1(ii)==ii || rndSeq1(ii)==rndBase(ii)
            rndSeq1(ii) = ceil(rand()*psExt);
        end
        while rndSeq2(ii)==rndSeq1(ii) || rndSeq2(ii)==rndBase(ii) || rndSeq2(ii)==ii
            rndSeq2(ii) = ceil(rand()*psExtB);
        end
    end
    
    posx = [pos;A];
    %     posx1 = [pos;pmid2];
    posr = pos(rndBase,:);
    posx1 = [posr;B];
    posxr = posx(rndSeq1,:);
    posxr1 = posx1(rndSeq2,:);
    %     posrr = pos(rndBase1,:);
    
    %%
    memory_rand_index1 = ceil(memory_size*rand(ps,1));
    %         memory_rand_index2 = ceil(memory_size*rand(ps,1));
    MUF = memory_MUF(memory_rand_index1);
    MUCr = memory_MUCr(memory_rand_index1);
    %%for generating crossover rate Cr
    %         Cr = normrnd(MUCr,0.1);
    Cr = laplacernd(MUCr,0.1,ps,1);
    term_Cr = find(MUCr == -1);
    Cr(term_Cr) = 0;
    Cr = min(Cr,1);
    Cr = max(Cr,0);
    
    label=zeros(ps,D);
    rndVal = rand(ps,D);
    onemat = zeros(ps,D);
    for ii = 1:ps
        label(ii,:) = rndVal(ii,:)<=Cr(ii);
        indexJ = ceil(rand()*D);
        onemat(ii,indexJ) = 1;
    end
    label = label|onemat;
    
    %% for generating scal factor F
    %     F = randCauchy(MUF,0.1);
    %     F = cos(g*pi/2*D)*MUF;
    % F = laplacernd(MUF,0.2,ps,1);
    % if g==1
    %     F = randCauchy(MUF,0.1);
    % else
    F = MUF;
    maxF = max(F);
    minF = min(F);
    if maxF~=minF
        for i=1:ps
            Fd(i)=(F(i)-minF)/(maxF-minF);
        end
        %         numf = numel(find((Fd<0.1)==1));
        F(Fd<0.1) = randCauchy(MUF(Fd<0.1),0.1);
        F(Fd>=0.1) = MUF(Fd>=0.1);
    else
        F = randCauchy(MUF,0.1);
    end
    % end
    
    term_F = find(F <= 0);
    while ~ isempty(term_F)
        F(term_F) = randCauchy(MUF(term_F),0.1);
        % F(term_F) = laplacernd(MUF(term_F),0.2,numel(term_F),1);
        term_F = find(F <= 0);
    end
    F = min(F, 1);
    
    
    F = max(F,0);
    gp = (max_nfes-nfes+1)/max_nfes;
    %     gp1 = (nfes)/max_nfes;
    %     w1=log10(10-gp1);
    w2=log(2-gp);
    pos = pbest + F.*(pbestB-pbest)+F.*(posxr1-posxr);
    %     pos = pbest + F.*(pbestB-pbest)+F.*(posr-posxr);
    pos(~label) = pbest(~label);
    %%
    
    
    pos = ((pos>=Xmin)&(pos<=Xmax)).*pos...
        +(pos<Xmin).*((Xmin+pbest).*rand(ps,D)/2) ...
        +(pos>Xmax).*((Xmax+pbest).*rand(ps,D)/2);
    
    dis = (pos-pbest).*label;
    
    
    eval = feval(fhd,pos',varargin{:});
    nfes = nfes + ps;
    bin = (pbestval > eval)';
    
    
    
    A = [A;pbest(bin==1,:)];
    lengthAdd = numel(pbest(bin==1));
    decayA = [decayA;zeros(lengthAdd,1)];
    if g==1
        ininum=numel(find(bin==1));
    end
    num_bin = numel(find(bin==1));
    decayr=num_bin/ininum;
    decayA = decayA + decayr;
    %     decayA = decayA +decayRate;
    %     if decayr == 0
    %     dcount = dcount + 1;
    %     end
    %                 if(numel(decayA)>0)
    %                     maxDecayA = max(decayA);
    %                 else
    %                     maxDecayA = 0;
    %                 end
    
    if size(A,1)>round(Archfactor*ps) ||decayr==0
        MergeA = [A,decayA];
        indexDecay = (decayA>1);
        MergeA(indexDecay,:) = [];
        len = length(MergeA(:,1));
        if len>round(Archfactor*ps)
            rndSel = randperm(len)';
            rndSel = rndSel(round(Archfactor*ps)+1:len);
            MergeA(rndSel,:) = [];
        end
        A = MergeA(:,1:D);
        decayA = MergeA(:,D+1);
        %         sizeA = round(Archfactor*ps);
        %         if size(A,1)>sizeA
        %             rndSel = randperm(size(A,1));
        %             rndSel = rndSel(sizeA+1:size(A,1));
        %             A(rndSel,:)=[];
    end
    
    pbest(bin==1,:) = pos(bin==1,:);
    pbestval(bin==1) = eval(bin==1);
    pmid = pos(bin==0,:);
    pmidval = eval(bin==0);
    %
    meval = mean(pbestval);
    bin1 = (pmidval < meval)';
    
    
    if size(B,1)>round(ArchfactorB*ps)
        B = B(1:round(ArchfactorB*ps),:);
    end
    B = [B;pmid(bin1,:)];
    
    
    sizeB = round(ArchfactorB*size(B,1));
    if size(B,1)>sizeB
        rndSel = randperm(size(B,1));
        rndSel = rndSel(sizeB+1:size(B,1));
        B(rndSel,:)=[];
    end
    SuccF = F(bin==1);
    SuccCr = Cr(bin==1);
    dis = dis(bin==1,:);
    dis = std(dis')';
    dis = dis/sum(dis);
    num_Succ = numel(SuccCr);
    memory_Ind = zeros(ps,memory_size);
    memory_Succ = zeros(ps,memory_size);
    for ii = 1:memory_size
        for j = 1:ps
            if memory_rand_index1(j) == ii
                memory_Ind(j,ii) = 1;
            end
        end
        memory_Succ(:,ii) = bin .* memory_Ind(:,ii);
    end
    memory_Ind_num = sum(memory_Ind);
    memory_Succ_num = sum(memory_Succ);
    if num_Succ > 0
        for ii = 1:memory_size
            if (memory_Ind_num(ii) == 0)
                probStra(ii) = 0.01;
            else
                probStra(ii) = memory_Succ_num(ii)^2/((sum(memory_Succ_num)+sum(memory_Ind_num))*sum(memory_Succ_num))+0.01;
            end
        end
        
        loc = find(probStra == min(probStra));
        indexLoc = loc(max(round(rand()*numel(loc)),1));
        memory_MUF(indexLoc) = (dis'*(SuccF.^2))/(dis'*SuccF);
        
        if max(SuccCr) == 0 || memory_MUCr(memory_order) == -1
            memory_MUCr(memory_order) = -1;
        else
            memory_MUCr(memory_order) = (dis'*(SuccCr.^2))/(dis'*SuccCr);
        end
        memory_order = memory_order + 1;
        if memory_order > memory_size
            memory_order = 1;
        end
    end
    
    vpop = 1;
    for ii=1:D
        iqr =abs(max(pos(:,ii))-min(pos(:,ii)));
        vpop = sqrt(vpop*(iqr));
    end
    RDI = sqrt(vpop/vlim);
    
    counter(bin==0) = counter(bin==0) + 1;
    counter(bin~=0) = 0;
    
    [~,Index_sel] = sort(pbestval);
    if  RDI > 0.001
        for i = 1:ps
            if counter(i)>=n && i~=gbestid
                sort_num = find(i==Index_sel);
                sort_num_s = sort_num(1);
                pbest_s_num = max(ceil(rand()*(sort_num_s-1)),1);
                pbest_s = pbest(Index_sel(pbest_s_num),:);
                %                 ggbest(i,:)=repmat((1-ff),1,D).*popmbest(i,:)+repmat(ff,1,D).*nbest(i,:);
                %                 F_con = repmat(randCauchy(0.5,0.1),1,D);
                %                 F_con = exp(1)/sqrt(5).*pi^(-1/8).*(randn.^2).*exp(-randn.^2/2);
                %                 F_con = cos(g*pi/2*D)*0.5;
                % F_con = exp(1)/sqrt(5).*pi^(-1/8).*(randn.^2).*exp(-rand.^2/2);
                F_con= randCauchy(0.5,0.1);
                %                  while F_con<0
                % %                     F_con = exp(1)/sqrt(5).*pi^(-1/8).*(randn.^2).*exp(-rand.^2/2);
                % F_con= laplacernd(0.5,0.2,1,1);
                %                  end
                %                 pbestX= pbest(i,:) + F_con.*(pbest_s - pbest(i,:)+ggbest(i,:)-pbest(i,:));
                pbestX= pbest(i,:) + F_con.*(pbest_s-pbest(i,:));
                pbestvalX = feval(fhd,pbestX',varargin{:});
                nfes=nfes+1;
                if pbestvalX < pbestval(i)
                    pbest(i,:) = pbestX;
                    pbestval(i) = pbestvalX;
                    counter(i)=0;
                end
            end
        end
    else
        %           if numel(find(counter>=n*ones(ps,1)))>=ps/ps_WW
        
        r_index = randi([1 ps],1,ps_WW);
        New_Point = [];
        FitVector = [];
        for i = 1 : ps_WW
            [NewP, fit] = WW_Process(popLS(i,:),fid,g,BestPoint,lu,fhd);
            New_Point = [New_Point;NewP];
            FitVector = [FitVector,fit];
            nfes = nfes + 1;
        end
        
        for i = 1 : ps_WW
            if FitVector(i) < pbestval(r_index(i))
                pbestval (r_index(i)) = FitVector(i);
                pbest(r_index(i),:) = New_Point(i,:);
                counter(r_index(i))=0;
            end
            
        end
    end
    if nfes< 0.5*max_nfes
        plan_ps = ceil((ps_min - ps_ini) / (2/3 * max_nfes - ps_ini)^2 * (nfes - ps_ini)^2 + ps_ini);
    else
        plan_ps = floor((ps_min-1/3*ps_ini)/((1/3)*max_nfes)*(nfes-max_nfes)+ps_min);
    end
    
    if ps > plan_ps
        reduc_nums = ps - plan_ps;
        if ps - reduc_nums < ps_min
            reduc_nums = ps - ps_min;
        end
        
        ps = ps - reduc_nums;
    end
    
    [gbestval,gbestid] = min(pbestval);
    gbest = pbest(gbestid,:);
    
    
    pbest(indexSel(ps+1:length(indexSel))',:) = [];
    counter(indexSel(ps+1:length(indexSel))',:) = [];
    pos = pbest;
    pbestval(indexSel(ps+1:length(indexSel))') = [];
    F(indexSel(ps+1:length(indexSel))') = [];
    [~,indexSel] = sort(pbestval);
    
    pbest_rate = pbest_rate_max+(pbest_rate_min-pbest_rate_max)*nfes/max_nfes;
    %      ArchfactorB = ArchfactorB_max-(ArchfactorB_max-ArchfactorB_min)*nfes/max_nfes;
    g = g +1;
    %          if((mod(nfes-ps,100*5*D)>=100*5*D-ps)&& mod(nfes,100*5*D)<ps)
    %         if mod(gen,5*D)==0
    if nfes>=print_num_record
        fprintf(fout,'%d\t%.15f\n',nfes,gbestval-targetbest(fid));
        print_num_record = print_num_record + interval_num;
    end
end
RecordT = toc;
if WriteFlag
    fclose(fout);
end
end
function result = randCauchy(mu, sigma)
[m,n] = size(mu);
result = mu + sigma*tan(pi*(rand(m,n)-0.5));
end
function x = laplacernd(a,b,m,n)
% Reference: http://en.wikipedia.org/wiki/Laplace_distribution
u = rand(m,n)-0.5;
x = a - b*sign(u).*log(1-2*abs(u));
end