%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 本程序为综合粒子群函数，其完全脱胎于mathworks的particleswarm函数，仅仅结构上与particleswarm类似，设置等操作不需要依赖于particleswarm
% 
% 暂时不支持PlotFcn绘图，目前只支持将每一次迭代中的最优解单独提出自行plot绘图
% 由于没有PlotFcn，本函数通过调节"参数设置.xlsx"的'基础部分'中Foreclosure进行强制停止行为，当Foreclosure在任意一次迭代中由FALSE转变成TRUE时，粒子群将在该迭代结束后终止
% 考虑到可能会设置一个外层嵌套for循环来多次运行粒子群算法，特加入一个保险迭代(Insurance
% iteration)来保证避免出现来不及切换Foreclosure而导致直接退出情况，可以在"参数设置.xlsx"中进行修改
% 
% 综合粒子群算法具有几种在particleswarm函数上并不存在的优化方式，分别为：交叉操作(Cross)、模拟退火接受函数(Anneal)、叠加态变异(Superposition)以及信赖域变异(Trust_region)
% 以上几种优化方式的参数都可以在"参数设置.xlsx"的Sheet中进行调整
% 如果需要关闭优化方式，可以去"参数设置.xlsx"中的'基础部分'的Judge进行调整，当Judge为FALSE时，将不会启用该优化。
% 几种优化方式中：
% 交叉为对于速度的变异
% 接受函数的温度并非模拟退火温度，而是一种与迭代次数相挂钩的温度，具体表示为Max_Iteration - Iteration + k，其中k为一个常数，可以在"参数设置"进行调整
% 叠加态变异是一种和云模型理论相类似的变异方式，具体原理可以参考：论正态云模型的普适性_李德毅
% 信赖域变异部分偷懒使用fmincon函数进行代笔
% 两种变异操作策略均为历史最优个体必然进行一次变异，剩下的变异位置留给随机粒子，并且会对速度进行一定的改变，具体变动见函数
% 
% 每种优化方式都配备了开始迭代，其中交叉和模拟退火接受函数还额外配备了结束迭代，可以在Sheet中进行调整
%
% 可以自己添加或者改进任意一优化方式，建议在初始化部分加入自己的新优化方式参数以便检查。
%
% 信赖域与叠加态两种优化方式由于是一种必良性的变异操作，需要在变异中计算结果，故每一次变异都具有很大的运算量，建议谨慎使用
%
% 本程序的运行时间与particleswarm几乎相同，需要开启并行时可以在"参数设置"的UseParallel中选择TRUE，另外信赖域变异也有并行选项，并行可以极大的节省代码运行时间
% 
% 本程序还额外支持以迭代次数作为除粒子所代表维度外的额外维进行动态演化的过程，如果需要开启则在Problemtype中选择TRUE，不开启时请务必关闭，若未关闭可能导致迭代次数占用一格输入向量位导致目标函数运行失败
% 
% 本程序对于原particleswarm的惯性权重方面进行了改动，现在每一个粒子都有对应的惯性，惯性更新公式见updateInertia函数
% 
% 本程序内并未内置混合函数，可以在函数外嵌套混合函数保证进入局部最优解
%
% Author：Ou Liu
% Email：202109031032@email.sxu.edu.cn 
% Phone number：(+86) 16670133062
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [x,fval,Best] = Syntheticpso(fun, nvars, lb, ub)
%% 初始化部分

  %% 初始粒子设置
InitialSwarmMatrix = readmatrix("参数设置",'Sheet','初始粒子');% 没有初始粒子会自动生成，生成方法为creationfcn
  %% options结构体设置
Syntheticpso = readmatrix('参数设置.xlsx','Sheet','基础部分');%读取粒子群参数
  %% 基础部分
Max_Iteration = Syntheticpso(1);%最大迭代次数
Min_Inertia = Syntheticpso(2);  %自适应惯性下界
Max_Inertia = Syntheticpso(3);  %自适应惯性上界
cSocial = Syntheticpso(4);      %社会学习因子
cSelf = Syntheticpso(5);        %自我学习因子
Swarsize = Syntheticpso(6);     %粒子个数
ObjectTarget = Syntheticpso(7); %理想函数值
  %% 优化部分
Judge_trust = Syntheticpso(8);          %是否启用信赖域变异
Judge_cross = Syntheticpso(9);          %是否启用交叉
Judge_anneal = Syntheticpso(10);        %是否启用退火接受函数
Judge_superposition = Syntheticpso(11);%是否启用叠加态变异
Min_neighbor = Syntheticpso(12);        %自适应邻域下界
    %% 信赖域变异设置
    trust_region = readmatrix('参数设置.xlsx','Sheet','trust_region');%读取信赖域变异参数

    MaxIter = trust_region(1);              %设置信赖域迭代次数
    Pm = trust_region(2);                   %设置变异概率
    Mut_Start_Iteration = trust_region(3);  %设置变异开始时间
    Tolcon = trust_region(4);               %设置信赖域判定epsilon收敛范围
    UseParallel = trust_region(5);          %设置是否启用并行

    trust_region = struct( ...                        设置信赖域变异结构体
        'Judge_trust_region',Judge_trust, ...         设置信赖域启用判定
        'MaxIteration',MaxIter, ...                   设置信赖域最大迭代次数
        'Pm',Pm, ...                                  设置信赖域变异概率
        'Mut_Start_Iteration',Mut_Start_Iteration, ...设置变异开始迭代
        'Tolcon',Tolcon, ...                           设置信赖域收敛判定范围
        'UseParallel',UseParallel ...                  设置是否启用并行
     );
    %% 交叉设置
    Cross = readmatrix('参数设置.xlsx','Sheet','Cross');%读取交叉参数

    Pc = Cross(1);             %设置交叉概率
    Start_Iteration = Cross(2);%设置交叉起始迭代
    End_Iteration = Cross(3);  %设置交叉终止迭代

    Cross = struct( ...                       设置交叉结构体
        'Judge_Cross',Judge_cross, ...        设置交叉启用判定
        'Pc',Pc, ...                          设置交叉率
        'Start_Iteration',Start_Iteration, ...设置交叉起始迭代
        'End_Iteration',End_Iteration ...     设置交叉终止迭代
        );
    %% 退火接受函数设置
    Anneal = readmatrix('参数设置.xlsx','Sheet','Anneal');%读取退火接受函数参数

    k = Anneal(1);              %设置模拟退火参数   
    Start_Iteration = Anneal(2);%设置接受起始迭代
    End_Iteration = Anneal(3);  %设置接受终止迭代

    Anneal = struct( ...                      设置接受函数结构体
        'Judge_Anneal',Judge_anneal, ...      设置接受函数启用判定
        'k',k, ...                            设置退火参数
        'Start_Iteration',Start_Iteration, ...设置接受起始迭代
        'End_Iteration',End_Iteration ...     设置接受终止迭代
        );
    %% 叠加态变异设置
    Superposition = readmatrix('参数设置.xlsx','Sheet',' Superposition');%读取叠加态变异参数

    numSuperposition = Superposition(1);  %设置叠加态个数
    sigma = Superposition(2);              %设置叠加态方差
    Pm = Superposition(3);                 %设置变异率
    Mut_Start_Iteration = Superposition(4);%设置变异开始迭代

    Superposition = struct( ...                        设置接受函数结构体
        'Judge_Superposition',Judge_superposition, ... 设置叠加态启用判定
        'numSuperposition',numSuperposition, ...       设置叠加态个数
        'sigma',sigma, ...                             设置叠加态方差
        'Pm',Pm, ...                                   设置变异率
        'Mut_Start_Iteration',Mut_Start_Iteration ...  设置变异起始迭代
        );
  %% 其他设置
UseParallel = Syntheticpso(13);         %是否启用并行
Problemtype = Syntheticpso(14);         %函数问题类型，如果为动态问题则为TRUE，静态问题则选择FALSE
Display = Syntheticpso(15);             %设置显示结果，如果为0则不显示，如果为1则显示最后一代，如果为2则显示若干代
InitialSwarmSpan = Syntheticpso(16);    %设置初始粒子范围
DisplayInterval = Syntheticpso(17);     %迭代次数间隔
StallIterLimit = Syntheticpso(18);      %显示窗口设置有关的参数(不是Plot)
Insurance_Iteration = Syntheticpso(19); %保险迭代窗口期设置


options = struct( ...% 设置结构体options
    'Max_Iteration',Max_Iteration, ...              设置最大迭代次数Max_Iteration
    'InertiaRange',[Min_Inertia,Max_Inertia], ...   设置自适应惯性上下界
    'cSocial',cSocial, ...                          设置社会学习因子
    'cSelf',cSelf, ...                              设置自我学习因子
    'Swarmsize',Swarsize, ...                       设置粒子个数
    'initialSwarmSpan',InitialSwarmSpan,...         设置粒子范围
    'initialSwarmMatrix',InitialSwarmMatrix,...     设置初始粒子
    'StallIterLimit',StallIterLimit,...             设置显示窗口设置有关的参数(不是Plot)
    'ObjectTarget',ObjectTarget, ...                设置目标函数值
    'Trust_region',trust_region, ...                设置信赖域变异
    'Cross',Cross, ...                              设置交叉
    'Anneal',Anneal,...                             设置接受函数
    'Superposition', Superposition,...              设置叠加态
    'Min_neighbor',Min_neighbor, ...                设置最小邻域半径
    'UseParallel',UseParallel, ...                  设置是否启用并行
    'Problemtype',Problemtype, ...                  设置问题类型
    'Verbosity',Display, ...                        设置显示格式    
    'DisplayInterval',DisplayInterval, ...          设置输出间隔
    'Insurance_Iteration',Insurance_Iteration ...   设置保险迭代期
    );
[x,fval,Best] = Core(fun,nvars,lb,ub,options);
end


function [x,fval,Best] = Core(fun,nvars,lb,ub,options)
%% 赋值阶段
pIdx = 1:options.Swarmsize;
exitFlag=[];
numParticles = options.Swarmsize;
minInertia = options.InertiaRange(1);
maxInertia = options.InertiaRange(2);
minNeighborhoodSize = max(2,floor(numParticles*options.Min_neighbor));
adaptiveNeighborhoodSize = minNeighborhoodSize;
%% 初始化阶段
state = makeState(nvars,lb,ub,fun,options);
bestFvals = min(state.Fvals);
bestFvalsWindow = nan(options.StallIterLimit, 1);

adaptiveInertiaCounter = 0;
%% 设置初始化惯性
if all(options.InertiaRange >= 0)
    adaptiveInertia = maxInertia*ones(numParticles,1);
elseif all(options.InertiaRange <= 0)
    adaptiveInertia = minInertia*ones(numParticles,1);
else
end

%% 设置显示部分
% Setup display header
if  options.Verbosity > 1
    fprintf('\n                                 Best            Mean     Stall\n');
    fprintf(  'Iteration     f-count            f(x)            f(x)    Iterations\n');
    fprintf('%5.0f         %7.0f    %12.4g    %12.4g    %5.0f\n', ...
        0, state.FunEval, bestFvals, mean(state.Fvals), 0);
end


i = 1;%记录迭代次数
while isempty(exitFlag)
    state.Iteration = state.Iteration + 1;
        % 更新邻域向量
        bestNeighborIndex = generateBestNeighborIndex(state,adaptiveNeighborhoodSize,numParticles);
        
        % 更新速度与位置
        state = updateParticles(state,adaptiveInertia,bestNeighborIndex,pIdx,numParticles,nvars,lb,ub,fun,options);
        

        % 求解粒子
        if options.Problemtype == false% 判断是否为动态问题
            state.Fvals = slove(state.Positions, fun, 1, options);
        else
            state.Fvals = slove([state.Positions,state.Iteration*ones(numParticles,1)],fun,1,options);
        end

        % 更新
        state = updateState(state,numParticles,pIdx,options);
        
        bestFvalsWindow(1+mod(state.Iteration-1,options.StallIterLimit)) = min(state.IndividualBestFvals);
        % 自适应惯性更新
        [state,adaptiveInertiaCounter,bestFvals,adaptiveNeighborhoodSize,adaptiveInertia] = updateInertia(state,options,bestFvals,adaptiveInertiaCounter,adaptiveNeighborhoodSize,adaptiveInertia,numParticles,minNeighborhoodSize);
        
        % 结束判断
        [exitFlag] = stopCore(options,state,bestFvalsWindow);
        [Bestfval(i),Posi] = min(state.Fvals);
        Best(i) = min(state.IndividualBestFvals);
        BestPosition(i,:) =  min(state.Positions(Posi));
        i = i+1;

end
if options.Problemtype == false %如果为静态问题，则为个体历史最优
        [fval,indexBestFval] = min(state.IndividualBestFvals);
        x = state.IndividualBestPositions(indexBestFval,:);
else%如果为动态问题，则需要将每一代都提出来
        x = BestPosition;
        fval = Bestfval;
        Best = Bestfval;

end
end





%% 设置粒子更新函数
function state =  updateParticles(state,adaptiveInertia,bestNeighborIndex,pIdx,numParticles,nvars,lb,ub,fun,options)
% 速度更新
state.Velocities(pIdx,:) = updateVelocities(state,adaptiveInertia,bestNeighborIndex,pIdx,nvars,options);
% 位置更新
[state.Positions(pIdx,:),tfInvalid] = updatePositions(state,lb,ub,pIdx,numParticles,nvars,fun,options);
% 速度修正
if any(tfInvalid(:))
    state.Velocities(tfInvalid) = 0;
end
end
%% 设置速度更新函数
function newVelocities = updateVelocities(state,adaptiveInertia,bestNeighborIndex,pIdx,nvars,options)

  %% 生成随机因子
randSelf = rand(numel(pIdx),nvars);
randSocial = rand(numel(pIdx),nvars);
%% 读取旧速度
oldVelocities = state.Velocities(pIdx,:);
  %%  速度更新公式
newVelocities = adaptiveInertia.*oldVelocities + ...旧速度部分
    options.cSelf*randSelf.*(state.IndividualBestPositions(pIdx,:)-state.Positions(pIdx,:)) + ...自我认知部分
    options.cSocial*randSocial.*(state.IndividualBestPositions(bestNeighborIndex(pIdx), :)-state.Positions(pIdx,:));%社会认知部分

  %%   交叉阶段
  if options.Cross.Judge_Cross == true %判断是否进行交叉
       if state.Iteration >= options.Cross.Start_Iteration
          if state.Iteration <= options.Cross.End_Iteration
           newVelocities = Crossfunction(newVelocities,pIdx,nvars,options);%进行交叉操作
          end
       end
  end
  %% 忽视速度无穷的情况
tfInvalid = ~all(isfinite(newVelocities), 2);
newVelocities(tfInvalid) = oldVelocities(tfInvalid);
end
%% 设置交叉函数
function newVelocities = Crossfunction(newVelocities,pIdx,nvars,options)
CrossNumber = randperm(numel(pIdx),round(numel(pIdx)*options.Cross.Pc));%设置交叉个体
CrossNumber = sort(CrossNumber);%对其进行升序排序
for i = 1:numel(CrossNumber)-1
CrossPosition = randperm(nvars,2);%在粒子中选择两个点位对其中间部分进行交叉
CrossPosition = sort(CrossPosition);%对其进行升序排序
a = CrossPosition(1);
b = CrossPosition(2);
temp=newVelocities(CrossNumber(i),(a:b));
newVelocities(CrossNumber(i),(a:b))=newVelocities(CrossNumber(i+1),(a:b));
newVelocities(CrossNumber(i+1),(a:b))=temp;
end
end



%% 设置位置更新函数
function [newPositions,tfInvalid] = updatePositions(state,lb,ub,pIdx,numParticles,nvars,fun,options)
lbMatrix = repmat(lb,options.Swarmsize,1);
ubMatrix = repmat(ub,options.Swarmsize,1);
  %% 更新位置
  newPositions = state.Positions(pIdx,:) + state.Velocities(pIdx,:);
  %% 去除大小为无穷的位置
  tfInvalid = any(~isfinite(newPositions), 2);
  tfInvalidFull = false(numParticles, 1);
  tfInvalidFull(pIdx) = tfInvalid;
  newPositions(tfInvalid, :) = state.Positions(tfInvalidFull, :);

  %% 边界约束
  tfInvalidLB = newPositions < lbMatrix(pIdx,:);
  if any(tfInvalidLB(:))
    tfInvalidLBFull = false(numParticles,nvars);
    tfInvalidLBFull(pIdx,:) = tfInvalidLB;
    newPositions(tfInvalidLB) = lbMatrix(tfInvalidLBFull);
    
    tfInvalid = tfInvalidLBFull;
  else
    tfInvalid = false(numParticles,nvars);
  end

  tfInvalidUB = newPositions > ubMatrix(pIdx,:);
  if any(tfInvalidUB(:))
    tfInvalidUBFull = false(numParticles,nvars);
    tfInvalidUBFull(pIdx,:) = tfInvalidUB;
    newPositions(tfInvalidUB) = ubMatrix(tfInvalidUBFull);
    
    tfInvalid = tfInvalid | tfInvalidUBFull;
  end
  %% 变异阶段
  % 信赖域变异
  if options.Trust_region.Judge_trust_region == true%判断是否进行信赖域变异
     newPositions = Mut_Trust_region(state,fun,newPositions,pIdx,lb,ub,options,nvars);%信赖域变异
  end
  % 叠加态变异
  if options.Superposition.Judge_Superposition == true%判断是否进行叠加态变异
      newPositions = Mut_Superposition(state,nvars,fun,newPositions,pIdx,lb,ub,options);%进行叠加态变异
  end
end
%% 设置信赖域变异函数
    %% 动态问题时间不变约束
function [c,ceq] = nonlcon(par)
load I

c = [];
ceq = [par(end)-I.state.Iteration];
end

  function newPositions = Mut_Trust_region(state,fun,newPositions,pIdx,lb,ub,options,nvars)
if state.Iteration <= options.Trust_region.Mut_Start_Iteration
    Pm = 0;
else
    Pm = options.Trust_region.Pm;
end

 MutationNumber = randperm(numel(pIdx),round(numel(pIdx)*Pm));%设置变异个体
 if options.Problemtype == false
for i=1:numel(MutationNumber)% 每次变异都选择一个最优剩下几个变异名额选择随机个体
    [~,MinFvalPosition] = min(state.IndividualBestFvals);
    if MutationNumber(i) == MinFvalPosition
        continue
    elseif i == numel(MutationNumber)
        MutationNumber(i) = MinFvalPosition;
    else
        continue
    end
end
 else
     for i=1:numel(MutationNumber)% 每次变异都选择一个最优剩下几个变异名额选择随机个体
    [~,MinFvalPosition] = min(state.Fvals);
    if MutationNumber(i) == MinFvalPosition
        continue
    elseif i == numel(MutationNumber)
        MutationNumber(i) = MinFvalPosition;
    else
        continue
    end
     end
end

 MutationNumber = sort(MutationNumber);%对其进行升序排序

     options.Trust_region.UseParallel = logical( options.Trust_region.UseParallel);

    opts = optimoptions('fmincon', ...%
    'TolCon',options.Trust_region.Tolcon, ...%
    'MaxFunEvals',options.Trust_region.MaxIteration, ...%
    'UseParallel', options.Trust_region.UseParallel, ...%
    'Algorithm','sqp', ...%
    'Display','none' ...%
    );
if options.Problemtype == false%静态无时变问题
    for i=1:numel(MutationNumber)% 开始信赖域自搜索变异
       [newPositions(MutationNumber(i),:),~,~,~,~,Velocities,~] = fmincon( fun,newPositions(MutationNumber(i),:),[],[],[],[],lb,ub,[],opts);
       newVelocities =Velocities';
       state.Velocities(MutationNumber(i),:) = newVelocities;
    end
else%动态问题，加入state.Iteration作为时间变量
     I = struct('nvars',nvars,'state',state);
     save I
    for i=1:numel(MutationNumber)% 开始信赖域自搜索变异
        [Positions,~,~,~,~,Velocities,~] = fmincon( fun,[newPositions(MutationNumber(i),:),state.Iteration],[],[],[],[],[lb,0],[ub,options.Max_Iteration],@nonlcon,opts);
        newPositions(MutationNumber(i),:) = Positions(1:end-1);
        Velocities = Velocities';
        state.Velocities(MutationNumber(i),:) = Velocities(1:end-1);

    end
end
end


%% 设置叠加态变异函数
function newPositions = Mut_Superposition(state,nvars,fun,newPositions,pIdx,lb,ub,options)
lbMatrix = repmat(lb,options.Swarmsize,1);
ubMatrix = repmat(ub,options.Swarmsize,1);
if state.Iteration <= options.Superposition.Mut_Start_Iteration
    Pm = 0;
else
    Pm = options.Superposition.Pm;
end


sigma = options.Superposition.sigma;
    if numel(sigma) >=2
    else
        sigma = sigma*(ones(1,nvars));
    end


numSuperposition = options.Superposition.numSuperposition;


 MutationNumber = randperm(numel(pIdx),round(numel(pIdx)*Pm));%设置变异个体


 if options.Problemtype == false
for i=1:numel(MutationNumber)% 每次变异都选择一个最优剩下几个变异名额选择随机个体
    [~,MinFvalPosition] = min(state.IndividualBestFvals);
    if MutationNumber(i) == MinFvalPosition
        continue
    elseif i == numel(MutationNumber)
        MutationNumber(i) = MinFvalPosition;
    else
        continue
    end
end
 else
     for i=1:numel(MutationNumber)% 每次变异都选择一个最优剩下几个变异名额选择随机个体
    [~,MinFvalPosition] = min(state.Fvals);
    if MutationNumber(i) == MinFvalPosition
        continue
    elseif i == numel(MutationNumber)
        MutationNumber(i) = MinFvalPosition;
    else
        continue
    end
     end
end
   for i=1:numel(MutationNumber)% 开始叠加态变异
       if options.UseParallel == true
       for j=1:numSuperposition
           if j == 1
               Superposition(j,:) = newPositions(MutationNumber(i),:);% 保留原位置
          else
               Superposition(j,:) = normrnd(newPositions(MutationNumber(i),:) ,sigma,1,nvars);% 将其他位置以sigma作为方差升级为叠加态位置
           end
       end
       else
           for j=1:numSuperposition
           if j == 1
               Superposition(j,:) = newPositions(MutationNumber(i),:);% 保留原位置
          else
               Superposition(j,:) = normrnd(newPositions(MutationNumber(i),:) ,sigma,1,nvars);% 将其他位置以sigma作为方差升级为叠加态位置
           end
           end
       end
       a = linspace(1,numSuperposition,numSuperposition); % 对于叠加态粒子个数进行拆分
       Invalid = any(~isfinite(Superposition), 2);
       tfInvalidFull = false(numSuperposition, 1);
       tfInvalidFull(a) = Invalid;
       Superposition(Invalid, :) = state.Positions(tfInvalidFull, :);
        lbMatrix = repmat(lbMatrix(1,:),numSuperposition,1);
        ubMatrix = repmat(ubMatrix(1,:),numSuperposition,1);
% Enforce bounds on positions and return logical array to update velocities where position exceeds bounds.
       tfInvalidLB = Superposition < lbMatrix(a,:);
      if any(tfInvalidLB(:))
        tfInvalidLBFull = false(numSuperposition,nvars);
         tfInvalidLBFull(a,:) = tfInvalidLB;
        Superposition(tfInvalidLB) = lbMatrix(tfInvalidLBFull);
        Invalid = tfInvalidLBFull;
      else
        Invalid = false(numSuperposition,nvars);
      end
   tfInvalidUB = Superposition > ubMatrix(1,:);
  if any(tfInvalidUB(:))
    tfInvalidUBFull = false(numSuperposition,nvars);
    tfInvalidUBFull(a,:) = tfInvalidUB;
    Superposition(tfInvalidUB) = ubMatrix(tfInvalidUBFull);
  end


  if options.Problemtype ==false
      Fvals(a) = slove( Superposition, fun, 1, options);  % 求出所有叠加态的解
  else
      Fvals(a) = slove([Superposition,state.Interation*ones(numSuperposition,1)],fun,1,options);

  end
        [~,Posi]=min(Fvals);%选取最优解，即进行观测坍缩
        state.Velocities(MutationNumber(i),:) = newPositions(MutationNumber(i),:) - Superposition(Posi,:);
        newPositions(MutationNumber(i),:)=Superposition(Posi,:);%最优解导回原位置

   end
end


%% 设置接受函数
function tfImproved = AnnealAcceptfunction(state,options)
for i=1:numel(state.Fvals)
    if state.Fvals(i) < state.IndividualBestFvals(i)
        tfImproved(i) = 1;
    else
      delE = state.Fvals(i) - state.IndividualBestFvals(i);
      h = 1/(1+exp(delE/(options.Max_Iteration - state.Iteration + options.Anneal.k)));%接受函数
            if h > rand
        tfImproved(i) = 1;
    else
        tfImproved(i) = 0;
            end
    end
end
end




%% 设置邻域更新函数
function bestNeighborIndex = generateBestNeighborIndex(state,adaptiveNeighborhoodSize,numParticles)

neighborIndex = zeros(numParticles, adaptiveNeighborhoodSize);%建立邻居矩阵

neighborIndex(:,1) = 1:numParticles; % 第一个邻居是自己

for i = 1:numParticles
    %% 确定排除粒子本身的随机邻居,即(numParticles-1)粒子
    neighbors = randperm(numParticles-1, adaptiveNeighborhoodSize-1);
    %% 对于>=当前粒子索引的索引加上1
    iShift = neighbors >= i;
    neighbors(iShift) = neighbors(iShift) + 1;
    neighborIndex(i,2:end) = neighbors;
end
%% 寻找邻域中的最优点
[~, bestRowIndex] = min(state.IndividualBestFvals(neighborIndex), [], 2);


bestLinearIndex = (bestRowIndex.'-1).*numParticles + (1:numParticles);
bestNeighborIndex = neighborIndex(bestLinearIndex);
end




%% 惯性、邻域范围更新函数
function [state,adaptiveInertiaCounter,bestFvals,adaptiveNeighborhoodSize,adaptiveInertia] = updateInertia(state,options,bestFvals,adaptiveInertiaCounter,adaptiveNeighborhoodSize,adaptiveInertia,numParticles,minNeighborhoodSize)
newBest = min(state.IndividualBestFvals);
if isfinite(newBest) && newBest < bestFvals
    bestFvals = newBest;
    state.LastImprovement = state.Iteration;
    state.LastImprovementTime = toc(state.StartTime);
    adaptiveInertiaCounter = max(0, adaptiveInertiaCounter-1);
    adaptiveNeighborhoodSize = minNeighborhoodSize;
else
    adaptiveInertiaCounter = adaptiveInertiaCounter+1;
    adaptiveNeighborhoodSize = min(numParticles, adaptiveNeighborhoodSize+minNeighborhoodSize);
end

% if adaptiveInertiaCounter < 2
%     adaptiveInertia = max(options.InertiaRange(1), min(options.InertiaRange(2), 2*adaptiveInertia));
% elseif adaptiveInertiaCounter > 5
%     adaptiveInertia = max(options.InertiaRange(2), min(options.InertiaRange(2), 0.5*adaptiveInertia));
% end
%% 更新自适应惯性权重范围
if options.UseParallel == true
parfor i = 1:options.Swarmsize
         if state.Fvals(i)<= mean(state.Fvals)%如果函数适应度小于均值
            adaptiveInertia(i)=options.InertiaRange(1)+ ...
            (state.Fvals(i) - min(state.Fvals))*(options.InertiaRange(2)-options.InertiaRange(1))/(mean(state.Fvals) - min(state.Fvals));
        else% 函数适应度大于等于均值
            adaptiveInertia(i) = options.InertiaRange(2);
        end 
end
else %不开启并行
for i = 1:options.Swarmsize
         if state.Fvals(i)<= meanf(state.Fvals)
            adaptiveInertia(i)=options.InertiaRange(1)+ ...
            (state.Fvals(i) - min(state.Fvals))*(options.InertiaRange(2)-options.InertiaRange(1))/(meanf(state.Fvals) - min(state.Fvals));
        else
            adaptiveInertia(i) = options.InertiaRange(2);
        end 
end
end
end

%% 结构体state更新最优位置与结果
function state = updateState(state,numParticles,pIdx,options)

state.FunEval = state.FunEval + numel(pIdx);

%% 接受函数
tfImproved = false(numParticles,1);
  if options.Anneal.Judge_Anneal == true%判断是否启用模拟退火接受函数
      if state.Iteration >= options.Anneal.Start_Iteration
          if state.Iteration <= options.Anneal.End_Iteration
      
          tfImproved = AnnealAcceptfunction(state,options);%模拟退火接受函数
          tfImproved = tfImproved';
          tfImproved = logical(tfImproved);
          end
      end

  else
    tfImproved(pIdx) = state.Fvals(pIdx) < ...
    state.IndividualBestFvals(pIdx);%原接受函数
  end
  %% 接受新解
  state.IndividualBestFvals(tfImproved) = state.Fvals(tfImproved);%更新值
  state.IndividualBestPositions(tfImproved, :) = state.Positions(tfImproved, :);%更新位置
end



%% 设置求解适应度函数
function y = slove(pop,fun,numObj,options)
 popSize = size(pop,1);
 y = zeros(popSize,numObj);
 if options.UseParallel ==true% 判断是否使用并行
 parfor (i = 1:popSize)
  y(i,:) = feval(fun,(pop(i,:)));
 end
 else
  for i = 1:popSize
  y(i,:) = feval(fun,(pop(i,:)));
  end
 end
end


%% 设置初始化函数
function options = initialParticleCheck(options)
numInitPositions = size(options.InitialSwarm, 1);

if numInitPositions == 0
    return
end

end 



%% 生成初始粒子函数
function swarm = CreationFcn(problemStruct)

nvars = problemStruct.nvars;
options = problemStruct.options;

[lb,ub] = determinePositionInitBounds(problemStruct.lb, problemStruct.ub, ...
    options.initialSwarmSpan);

numParticles = options.Swarmsize;
numInitPositions = size(options.initialSwarmMatrix, 1);
numPositionsToCreate = numParticles - numInitPositions;

%% 创建初始化粒子
swarm = zeros(numParticles,nvars);

%% 手动输入粒子部分
if numInitPositions > 0
    swarm(1:numInitPositions,:) = options.InitialSwarm;
end

% Create remaining particles, randomly sampling within lb and ub
span = ub - lb;
swarm(numInitPositions+1:end,:) = repmat(lb,numPositionsToCreate,1) + ...
    repmat(span,numPositionsToCreate,1) .* rand(numPositionsToCreate,nvars);

end


%% 确定初始位置上下界
function [lb,ub] = determinePositionInitBounds(lb,ub,initialSwarmSpan)
% Update lb and ub using positionInitSpan, so that initial bounds are
% always finite
lbFinite = isfinite(lb);
ubFinite = isfinite(ub);
lbInf = ~lbFinite;
ubInf = ~ubFinite;

% If lb and ub are both finite, do not update the bounds.

% If lb & ub are both infinite, center the range around 0.
idx = lbInf & ubInf;
lb(idx) = -initialSwarmSpan(idx)/2;
ub(idx) = initialSwarmSpan(idx)/2;

% If only lb is finite, start the range at lb.
idx = lbFinite & ubInf;
ub(idx) = lb(idx) + initialSwarmSpan(idx);

% If only ub is finite, end the range at ub.
idx = lbInf & ubFinite;
lb(idx) = ub(idx) - initialSwarmSpan(idx);
end





%% 创建初始结构体
function state = makeState(nvars,lb,ub,fun,options)
% 创建初始上下界
lbMatrix = repmat(lb,options.Swarmsize,1);
ubMatrix = repmat(ub,options.Swarmsize,1);
% 创建初始结构体
state = struct;
state.Iteration = 0; % current generation counter
state.StartTime = tic; % tic identifier
state.StopFlag = false; % OutputFcns flag to end the optimization
state.LastImprovement = 1; % generation stall counter
state.LastImprovementTime = 0; % stall time counter
state.FunEval = 0;
numParticles = options.Swarmsize;

%% 初始化
    problemStruct = struct( ...
        'objective',fun, ...
        'lb',lb, ...
        'ub',ub, ...
        'nvars',nvars, ...
        'options',options ...
        );
    state.Positions = CreationFcn(problemStruct);

%% 边界约束
if any(any(state.Positions < lbMatrix)) || any(any(state.Positions > ubMatrix))
    state.Positions = max(lbMatrix, state.Positions);
    state.Positions = min(ubMatrix, state.Positions);
end

% 随机采样生成速度
vmax = min(ub-lb, options.initialSwarmSpan);
state.Velocities = repmat(-vmax,numParticles,1) + ...
    repmat(2*vmax,numParticles,1) .* rand(numParticles,nvars);

if options.Problemtype == false %如果不带时变
    fvals = slove(state.Positions(1:end,:),fun,1,options);
    state.Fvals = fvals;
else% 带时变
    fvals = slove([state.Positions(1:end,:),state.Iteration*ones(options.Swarmsize,1)],fun,1,options);%加入一个和迭代次数相关的时间作为变量的一部分
    state.Fvals =  fvals;  
end

state.FunEval = numParticles;

state.IndividualBestFvals = state.Fvals;
state.IndividualBestPositions = state.Positions;
end 


%% 设置粒子群停止条件
function [exitFlag] = stopCore(options,state,bestFvalsWindow)
iteration = state.Iteration;

if options.Problemtype == true
   bestFval = min(state.Fvals);
else
    iterationIndex = 1+mod(iteration-1,options.StallIterLimit);
    bestFval = bestFvalsWindow(iterationIndex);
end



if options.Verbosity > 1 && ...
        mod(iteration,options.DisplayInterval)==0 && ...
        iteration > 0
    FunEval  = state.FunEval;
    MeanFval = meanf(state.Fvals);
    StallGen = iteration  - state.LastImprovement;
    fprintf('%5.0f         %7.0f    %12.4g    %12.4g    %5.0f\n', ...
        iteration, FunEval, bestFval, MeanFval, StallGen);
end
Zne = readmatrix('参数设置.xlsx','Sheet','基础部分');
Foreclosure = Zne(20);%设置强制退出
reasonToStop = '';
exitFlag = [];
if state.Iteration >= options.Insurance_Iteration
    if state.Iteration >= options.Max_Iteration
        exitFlag = 0;
    elseif Foreclosure == true %强制退出
        exitFlag = 1;
    end
end

if ~isempty(reasonToStop) && options.Verbosity > 0
    return
end

% 重新打印
if options.Verbosity > 1 && rem(iteration,30*options.DisplayInterval)==0 && iteration > 0
    fprintf('\n                                 Best            Mean     Stall\n');
    fprintf('Iteration     f-count            f(x)            f(x)    Iterations\n');
end

end



%% 均值函数
function m = meanf(x)
tfValid = ~isnan(x);
n = sum(tfValid);
if n==0
    % 避免出现m/0的情况
    m = NaN;
else
    m = sum(x(tfValid)) ./ n;
end
end 

