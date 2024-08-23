%%--------------�����Ż��㷨----------------------%%
%% ���룺
%   pop:��Ⱥ����
%   dim:�������ǵ�ά��
%   ub:�����ϱ߽���Ϣ��ά��Ϊ[1,dim];
%   lb:�����±߽���Ϣ��ά��Ϊ[1,dim];
%   fobj:Ϊ��Ӧ�Ⱥ����ӿ�
%   maxIter: �㷨�����������������ڿ����㷨��ֹͣ��
%% �����
%   Best_Pos��Ϊ�����㷨�ҵ�������λ��
%   Best_fitness: ����λ�ö�Ӧ����Ӧ��ֵ
%   IterCure:  ���ڼ�¼ÿ�ε����������Ӧ�ȣ��������������Ƶ������ߡ�
function [Best_Pos,Best_fitness,IterCurve] = GWO2(pop,dim,ub,lb,fobj,maxIter)
   
    %% ����Alpha��Beta��Delta��
    Alpha_pos=zeros(1,dim);
    Alpha_score=inf; 

    Beta_pos=zeros(1,dim);
    Beta_score=inf; 

    Delta_pos=zeros(1,dim);
    Delta_score=inf; 
    %% ��ʼ����Ⱥλ��
    Positions = initialization2(pop,ub,lb,dim);
    %% ������Ӧ��ֵ
    fitness = zeros(1,pop);
    for i = 1:pop
       fitness(i) = fobj(Positions(i,:));
    end
    %% ����Ӧ�������ҵ�Alpha��Beta��Delta��
    %Ѱ����Ӧ����С��λ��
    [SortFitness,indexSort] = sort(fitness);
    %��¼��Ӧ��ֵ��λ��
    Alpha_pos = Positions(indexSort(1),:);
    Alpha_score = SortFitness(1);
    Beta_pos = Positions(indexSort(2),:);
    Beta_score = SortFitness(2);
    Delta_pos = Positions(indexSort(3),:);
    Delta_score = SortFitness(3);
    gBest = Alpha_pos;
    gBestFitness = Alpha_score;  
    %��ʼ����
    for t = 1:maxIter
        %����a��ֵ
       %a=2-t*((2)/maxIter);      
        a_min=0;
        a_max=2;
        alpha=0.3;
         a=(a_max-a_min).*exp(-(t./(maxIter.*alpha)).^2)+a_min;%aΪ��������  
        for i = 1:pop
            for j = 1:dim
              %% ����Alpha�Ǹ���λ��
              
              a11=log(fobj(Alpha_pos).^2);
              a22=log(fobj(Beta_pos).^2);
              a33=log(fobj(Delta_pos).^2);
              
               w_X1=a11/(a11+a22+a33);
                w_X2=a22/(a11+a22+a33);
                w_X3=a33/(a11+a22+a33);
               
              
                r1=rand(); %[0,1]�����
                r2=rand(); % [0,1]�����            
                A1=2*a*r1-a; % ����A1
                C1=2*r2; % ����C1
                D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % �������Alpha�ľ���
                X1=w_X1*Alpha_pos(j)-A1*D_alpha; %λ�ø���
              
              %% ����Beta�Ǹ���λ��
                r1=rand();%[0,1]�����
                r2=rand();%[0,1]�����       
                A2=2*a*r1-a;% ����A2
                C2=2*r2; % ����C2           
                D_beta=abs(C2*Beta_pos(j)-Positions(i,j));  % �������Beta�ľ���
                X2=w_X2*Beta_pos(j)-A2*D_beta; %λ�ø���
                
              %% ����Delta�Ǹ���λ��
                r1=rand();%[0,1]�����
                r2=rand();%[0,1]�����
                A3=2*a*r1-a; % ����A3
                C3=2*r2; % ����C3 
                D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % �������Delta�ľ���
                X3=w_X3*Delta_pos(j)-A3*D_delta; %λ�ø���
                %����λ��
                Positions(i,j)=(X1+X2+X3)/3;
            end
            %% �߽���
             Positions(i,:) = BoundaryCheck(Positions(i,:),ub,lb,dim);
        end
        %������Ӧ��ֵ
         for i = 1:pop
             fitness(i) = fobj(Positions(i,:));
             % ���� Alpha, Beta,  Delta��
            if  fitness(i)<Alpha_score 
                Alpha_score= fitness(i); % ����alpha��
                Alpha_pos=Positions(i,:);
            end
            if  fitness(i)>Alpha_score &&  fitness(i)<Beta_score 
                Beta_score= fitness(i); % ����beta��
                Beta_pos=Positions(i,:);
            end
            if  fitness(i)>Alpha_score &&  fitness(i)>Beta_score &&  fitness(i)<Delta_score 
                Delta_score= fitness(i); %����delta��
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