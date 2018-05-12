%% Algoritmo Genético para otimização topológica
% Inicialmente gera uma população aleatória de candidatos
% Calcula a fitness (função objetiva a ser maximizada) dos indivíduos da
%   população original
% Os melhores indivíduos se intercruzam, com um crossover dos cromosomos
% Os cromosomos, então, recebem mutações aleatórias
% Alguns melhores da iteração passada se mantém inalterados, para garantir
%   que não haja piora na fitness
% O ciclo se repete, até atingir um número limite de iterações
% O melhor resultado de todos os ciclos é então escolhido.
%
% No caso, este algoritmo resolve um problema de compliance.
% Pode ser adaptado para outros casos.
% Há algumas maneiras de se fazer a seleção e cruzamento, como intercompetição e
% recompensas. Aqui, se usa probabilidade proporcional à fitness.
% 
% Normalmente, para a codificação dos genes, usa-se uma concatenação de 
% código binário para representar cada característica, mas aqui, como o
% problema é sólido-vazio, pode-se usar uma concatenação direta da matriz
% sólido-vazio.
%
%%%%%%%%%%%%%%%%%%%%%%Preciso encontrar uma maneira de manter somente
%%%%%%%%%%%%%%%%%%%%%%volfrac porcento dos genes expressos!


%% INICIALIZAÇÃO DAS VARIÁVEIS DO PROBLEMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nelx=30; nely=15;
% vol = 1; change = 1;
ij = 0;
% figx=8; figy=5; %%Quantidade de colunas e linhas a serem usadas no display da população em cada geração
%% DADOS GENÉTICOS
rng('shuffle') 
itlim=50;
maxfitness(1:itlim)=0; %%Melhor fitness de cada iteração
%bestfitness(1:itlim)=0; %%Melhor fitness global
meanfitness(1:itlim)=0; %%Fitness mediana de cada iteração
averagefitness(1:itlim)=0; %%Fitness média de cada iteração
mutationrate=0.20; %%Porcentagem do genoma que pode ser alterado em cada iteração
popsize=40; %%Tamanho da população
%chromosome=zeros(nelx*nely,popsize); %%Pré-alocamento dos genes, a ser esvaziado posteriormente
deathrate=0.05; %%Quantos porcento da população morre em cada geração - o complementar para 1 de quantos se mantém
deathn=round((1-deathrate)*popsize);
rep_probability=zeros(1,popsize); %%Porcentagem de um cromosomo ser escolhido para reprodução
passeddowngenes=0.6; %%Porcentagem de genes a serem passados do melhor pai para a criança
%passeddowngenes=round(passeddowngenes*nely*nelx);
pairs=zeros(popsize,2); %%Array de relação de casais
%% DADOS DO MEIO MATERIAL
%objective=@(x,y) 20+x^2+y^2-10*(2*cos(2*pi*x)+2*cos(2*pi*y)); %%Rastrigin
objective=@(x,y) (x^2+y-11)^2+(x+y^2-7)^2; %Himmelblau
%objective=@(x,y) -(y+47)*sin(sqrt(abs(x/2+y+47)))-x*sin(sqrt(abs(x-(y+47))));
dim=nargin(objective);
fitness(1:popsize)=0;
lim=[-5,5; ... %%Limites das variáveis do "habitat", predefinidos pelo usuário
      -5,5];
radius=abs(lim(1,2)-lim(1,1))/5; %%Magnitude da mutação induzida
step=0.1; %%Intervalo entre cada ponto calculado do meio %%PELO AMOR DE DEUS **NÃO** COLOQUE BAIXO DEMAIS
interpoints1=lim(1,1):step:lim(1,2);
interpoints2=lim(2,1):step:lim(2,2);
[X,Y]=meshgrid(interpoints1,interpoints2);
%objmatrix=-(Y+47).*sin(sqrt(abs(X./2+Y+47)))-X.*sin(sqrt(abs(X-(Y+47)))); % Egg
objmatrix=(X.^2+Y-11).^2+(X+Y.^2-7).^2; % Himmelblau
subplot(1,2,1)
habitat=surf(X,Y,objmatrix);
habitat.EdgeColor = 'none';
subplot(1,2,2)
contourf(X,Y,objmatrix)
%% INÍCIO DO PROGRAMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INICIALIZAR A POPULAÇÃO INICIAL
population= rand(popsize,dim);
for i=1:dim
    population(:,i) = lim(i,1)+(lim(i,2)-lim(i,1))*population(:,i);
end
%population(:,:,5)=ones(nely,nelx);
while ij<itlim
    clf();
    ij = ij + 1;
    %% ANÁLISE DE FITNESS DE CADA GERAÇÃO
    %%Encoding dos cromossomos
%     for individual=1:popsize
%         chromosome(:,individual)=([num2str(population(:,individual));
%     end
    %%CALCULO DA FUNÇÃO FITNESS
    % -Tipo, só checa a função nos pontos de cada iteração l m a o
    for individual=1:popsize
        fitness(individual)=objective(population(individual,1),population(individual,2));
    end
    [fitness,ordem]=sort(fitness,'ascend');
    maxfitness(ij)=fitness(1);
    if ij==1
        bestfitness(ij)=maxfitness(ij);
    else if ij~=1
            if maxfitness(ij)<min(bestfitness)
                bestfitness(ij)=maxfitness(ij);
            else
                bestfitness(ij)=bestfitness(ij-1);
            end
        end
    end
    meanfitness(ij)=(fitness(ceil(popsize/2))+fitness(floor(popsize/2)))/2;
    averagefitness(ij)=sum(fitness)/popsize;
    %%DISPLAY DOS INDIVÍDUOS COM FITNESS
    fprintf('Generation:%3i Generation Fitness:%8.4f Best Fitness:%8.4f \n',(ij),maxfitness(ij),bestfitness(ij));
    subplot(1,2,1)
    habitat=surf(interpoints1,interpoints2,objmatrix,'FaceAlpha',0.5);
    habitat.EdgeColor = 'none';
    hold on
    for i=1:popsize
        scatter3(population(i,1),population(i,2),fitness(i),70,'filled');
    end
    title(['Habitat. Current best Fitness: ' num2str(maxfitness(ij)) ', mean fitness: ' num2str(meanfitness(ij))]);
    colormap bone;
    subplot(1,2,2)
    contourf(X,Y,objmatrix)
    colormap bone;
    hold on
    population=population(ordem,:);
    for i=1:popsize
        scatter(population(i,1),population(i,2),70,'filled');
    end
    colorbar
    title(['Top-down habitat. Overall best Fitness: ' num2str(bestfitness(ij)) ', avg fitness: ' num2str(averagefitness(ij))]);
    %chromosome=chromosome(:,ordem);
    pause(1)
    %% SELEÇÃO DOS INDIVÍDUOS
    %%CALCULO DE PROBABILIDADES DE REPRODUCAO
    
    %filter=sum(fitness); %%Fitness proportionate selection
    
    tempfitness=fitness;%%Fitness proportionate selection corrigido (3 linhas)
    tempfitness=tempfitness+abs(maxfitness(1));
    filter=sum(tempfitness);
    
    %filter=sum(1:popsize); %%Normalized fitness proportionate selection (max chance: popsize/sum(1:popsize) )
    for i=1:popsize
        %rep_probability(i)=fitness(i)/filter; %Fitness proportionate
        %rep_probability(i)=(popsize-i+1)/sum(1:popsize); %Normalized, shitty selection
        rep_probability(i)=1-(fitness(i)+50)/filter;  %% eu não sei porque isso funciona, mas funciona
        cumprob(i)=sum(rep_probability(1:i));
    end
    
    %%SELEÇÃO DOS PAIS
    for i=1:popsize
        lucky=rand;
        j=1;
         %%SELEÇÃO DO PRIMEIRO PARENTE
        if lucky>=cumprob(j)
            while lucky>=cumprob(j)
                j=j+1;
                if lucky<cumprob(j)
                    pairs(i,1)=j;
                    break
                end
            end
        else
            pairs(i,1)=j;
        end
        %%SELEÇÃO DO SEGUNDO PARENTE
        % NORMAL, SEM PREVENÇÃO DE INCESTO
        lucky=rand;
        j=1;
        if lucky>=cumprob(j)
            while lucky>=cumprob(j)
                j=j+1;
                if lucky<cumprob(j) %%exato mesmo código pq preguiça
                    pairs(i,2)=j;
                    break
                end
            end
        else
            pairs(i,2)=j;
        end
        
        % PREVENÇÃO DE INCESTO CONSIGO MESMO (EXPERIMENTAL)
        % selmethod=randi([0,rep_probability(j)]);
        % tempprob=cumprob-rep_probability(j);
        
        % CHECAGEM DE DOMINÂNCIA E CORREÇÃO
        if fitness(pairs(i,2))<fitness(pairs(i,1))
            temp=pairs(i,2);
            pairs(i,2)=pairs(i,1);
            pairs(i,1)=temp;
        end
    end
    %% CROSSOVER DOS GENES
    oldpopulation=population;
%     oldchromosome=chromosome;
    %%PAR RECESSIVO COMO BASE, SERÁ ACRESCENTADO OS GENES DOMINANTES
%     for i=1:popsize
%         crossover(:,i)=chromosome(:,pairs(i,2)); %%Código com encoding
%     end

%     for i=1:popsize
%         crossover(:,i)=fitness(pairs(i,2)); %%Sem encoding
%     end

%         dominantgenes=datasample(1:nely*nelx,passeddowngenes,'replace',false); %dos genes dominantes serem passados
%         passeddowngenes
%         dominantgenes
%         crossover(:,1)
%         figure(2)
%         colormap(gray); imagesc(-population(:,:,pairs(1,2))); axis equal; axis tight; axis off;
%         figure(3)
%         for j=1:passeddowngenes 
%             crossover(dominantgenes(j),1)=chromosome(dominantgenes(j),pairs(1,1));
%         end
%         test=reshape(crossover(:,1),[nely,nelx]);
%         
%         colormap(gray); imagesc(-test); axis equal; axis tight; axis off;
        
%         for j=1:passeddowngenes 
%             crossover(dominantgenes(j),1)=chromosome(dominantgenes(j),pairs(1,1));
%         end
%         figure(3)
%         olormap(gray); imagesc(-population(:,:,pairs(1,2))); axis equal; axis tight; axis off;
    %%CROSSOVER COM OS GENES DOMINANTES
%     for i=1:popsize
%         rng('shuffle')                                                         %Garantia de passeddowngenes%
%         dominantgenes=datasample(1:nely*nelx,passeddowngenes,'replace',false); %dos genes dominantes serem passados
%         for j=1:passeddowngenes 
%             crossover(dominantgenes(j),i)=chromosome(dominantgenes(j),pairs(i,1));
%         end
%         population(:,:,i)=reshape(crossover(:,i),[nely,nelx]);
%         chromosome(:,i)=crossover(:,i);
%     end
    for i=1:popsize                                                       
        for j=1:dim
            population(i,j)=passeddowngenes*oldpopulation(pairs(i,1),j)+(1-passeddowngenes)*oldpopulation(pairs(i,2),j);
        end
    end
    %% MUTATIONS
    %%DETECÇÃO SE ESTÁ CONVERGINDO DEMAIS
%     if ij>10
%         if abs(sum(maxfitness(ij-9:ij-5))-sum(maxfitness(ij-4:ij)))/sum(maxfitness(ij-4:ij))<0.3
%             mutation=mutation+0.05;
%         else
%             mutation=mutation-0.05;
%         end
%     end
%     %%MUTAÇÃO BITFLIP SIMPLES
%     for i=1:popsize
%         mutation=rand;
%         %%LEITURA DOS CROMOSOMOS
%         for m=1:(nelx*nely)
%             if mutation<mutationrate
%                 %%MUTAÇÃO
%                 if chromosome(m,i)==1
%                     chromosome(m,i)=0;
%                 else
%                     chromosome(m,i)=1;
%                 end
%             end
%         end
%         %%REFORMAÇÃO
%         population(:,:,i)=reshape(chromosome(:,i),[nely,nelx]);
%     end

%     %%MUTAÇÃO FLOATING POINT SIMPLES
    for i=1:popsize
        for j=1:dim
            if rand<mutationrate
                mutation=radius*2*(rand-0.5);
                population(i,j)=population(i,j)+mutation;
            end
        end
        
    end

    %% ELITISMO, parece não funcionar, pelo menos em funções simples (ex: 2d)
%     for sin=1:deathn
%         population(popsize-sin,:)=oldpopulation(sin,:);
%     end
end

figure()
%objmatrix=-(Y+47).*sin(sqrt(abs(X./2+Y+47)))-X.*sin(sqrt(abs(X-(Y+47)))); %egg

habitat=surf(interpoints1,interpoints2,objmatrix,'FaceAlpha',0.5);
habitat.EdgeColor = 'none';
colormap bone;
hold on
for i=1:popsize
    scatter3(oldpopulation(i,1),oldpopulation(i,2),fitness(i),70,'filled');
end
colorbar
size(maxfitness)
size(averagefitness)
size(meanfitness)
size(bestfitness)
title(['Final Habitat. Best Fitness: ' num2str(bestfitness(ij-1)) ', mean fitness: ' num2str(meanfitness(ij-1)) ', average fitness: ' num2str(averagefitness(ij-1))]);
pause(1e-6)
figure()
title(['Fitness statistics.']);
plot(1:ij,bestfitness,'-kd')
hold on
plot(1:ij,maxfitness,'-bo')
plot(1:ij,meanfitness,'-gs')
plot(1:ij,averagefitness,'-r^')
legend('Best overall fitness','Best generation fitness','Mean generation fitness','Average generation fitness')
xlabel('Generations')
ylabel('Fitness coefficients')
grid on
grid minor