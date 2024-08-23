%%--------------灰狼优化算法----------------------%%
%% 输入：
%   pop:种群数量
%   dim:单个灰狼的维度
%   ub:灰狼上边界信息，维度为[1,dim];
%   lb:灰狼下边界信息，维度为[1,dim];
%   fobj:为适应度函数接口
%   maxIter: 算法的最大迭代次数，用于控制算法的停止。
%% 输出：
%   Best_Pos：为灰狼算法找到的最优位置
%   Best_fitness: 最优位置对应的适应度值
%   IterCure:  用于记录每次迭代的最佳适应度，即后续用来绘制迭代曲线。
function [Best_Pos,Best_fitness,IterCurve] = GWO2(pop,dim,ub,lb,fobj,maxIter)
   
    %% 定义Alpha，Beta，Delta狼
    Alpha_pos=zeros(1,dim);
    Alpha_score=inf; 

    Beta_pos=zeros(1,dim);
    Beta_score=inf; 

    Delta_pos=zeros(1,dim);
    Delta_score=inf; 
    %% 初始化种群位置
    Positions = initialization2(pop,ub,lb,dim);
    %% 计算适应度值
    fitness = zeros(1,pop);
    for i = 1:pop
       fitness(i) = fobj(Positions(i,:));
    end
    %% 对适应度排序，找到Alpha，Beta，Delta狼
    %寻找适应度最小的位置
    [SortFitness,indexSort] = sort(fitness);
    %记录适应度值和位置
    Alpha_pos = Positions(indexSort(1),:);
    Alpha_score = SortFitness(1);
    Beta_pos = Positions(indexSort(2),:);
    Beta_score = SortFitness(2);
    Delta_pos = Positions(indexSort(3),:);
    Delta_score = SortFitness(3);
    gBest = Alpha_pos;
    gBestFitness = Alpha_score;  
    %开始迭代
    for t = 1:maxIter
        %计算a的值
       %a=2-t*((2)/maxIter);      
        a_min=0;
        a_max=2;
        alpha=0.3;
         a=(a_max-a_min).*exp(-(t./(maxIter.*alpha)).^2)+a_min;%a为收敛因子  
        for i = 1:pop
            for j = 1:dim
              %% 根据Alpha狼更新位置
              
              a11=log(fobj(Alpha_pos).^2);
              a22=log(fobj(Beta_pos).^2);
              a33=log(fobj(Delta_pos).^2);
              
               w_X1=a11/(a11+a22+a33);
                w_X2=a22/(a11+a22+a33);
                w_X3=a33/(a11+a22+a33);
               
              
                r1=rand(); %[0,1]随机数
                r2=rand(); % [0,1]随机数            
                A1=2*a*r1-a; % 计算A1
                C1=2*r2; % 计算C1
                D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % 计算距离Alpha的距离
                X1=w_X1*Alpha_pos(j)-A1*D_alpha; %位置更新
              
              %% 根据Beta狼更新位置
                r1=rand();%[0,1]随机数
                r2=rand();%[0,1]随机数       
                A2=2*a*r1-a;% 计算A2
                C2=2*r2; % 计算C2           
                D_beta=abs(C2*Beta_pos(j)-Positions(i,j));  % 计算距离Beta的距离
                X2=w_X2*Beta_pos(j)-A2*D_beta; %位置更新
                
              %% 根据Delta狼更新位置
                r1=rand();%[0,1]随机数
                r2=rand();%[0,1]随机数
                A3=2*a*r1-a; % 计算A3
                C3=2*r2; % 计算C3 
                D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % 计算距离Delta的距离
                X3=w_X3*Delta_pos(j)-A3*D_delta; %位置更新
                %更新位置
                Positions(i,j)=(X1+X2+X3)/3;
            end
            %% 边界检查
             Positions(i,:) = BoundaryCheck(Positions(i,:),ub,lb,dim);
        end
        %计算适应度值
         for i = 1:pop
             fitness(i) = fobj(Positions(i,:));
             % 更新 Alpha, Beta,  Delta狼
            if  fitness(i)<Alpha_score 
                Alpha_score= fitness(i); % 更新alpha狼
                Alpha_pos=Positions(i,:);
            end
            if  fitness(i)>Alpha_score &&  fitness(i)<Beta_score 
                Beta_score= fitness(i); % 更新beta狼
                Beta_pos=Positions(i,:);
            end
            if  fitness(i)>Alpha_score &&  fitness(i)>Beta_score &&  fitness(i)<Delta_score 
                Delta_score= fitness(i); %更新delta狼
                Delta_pos=Positions(i,:);
            end
         end
        gBest = Alpha_pos;
        gBestFitness = Alpha_score;  
        IterCurve(t) = gBestFitness;  
    end
    Best_Pos = gBest;
    Best_fitness = gBestFitness;
end