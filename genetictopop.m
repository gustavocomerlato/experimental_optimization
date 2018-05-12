%% Algoritmo Gen�tico para otimiza��o topol�gica
% Inicialmente gera uma popula��o aleat�ria de candidatos
% Calcula a fitness (fun��o objetiva a ser maximizada) dos indiv�duos da
%   popula��o original
% Os melhores indiv�duos se intercruzam, com um crossover dos cromosomos
% Os cromosomos, ent�o, recebem muta��es aleat�rias
% Alguns melhores da itera��o passada se mant�m inalterados, para garantir
%   que n�o haja piora na fitness
% O ciclo se repete, at� atingir um n�mero limite de itera��es
% O melhor resultado de todos os ciclos � ent�o escolhido.
%
% No caso, este algoritmo resolve um problema de compliance.
% Pode ser adaptado para outros casos.
% H� algumas maneiras de se fazer a sele��o e cruzamento, como intercompeti��o e
% recompensas. Aqui, se usa probabilidade proporcional � fitness.
% 
% Normalmente, para a codifica��o dos genes, usa-se uma concatena��o de 
% c�digo bin�rio para representar cada caracter�stica, mas aqui, como o
% problema � s�lido-vazio, pode-se usar uma concatena��o direta da matriz
% s�lido-vazio.
%
%%%%%%%%%%%%%%%%%%%%%%Preciso encontrar uma maneira de manter somente
%%%%%%%%%%%%%%%%%%%%%%volfrac porcento dos genes expressos!


%% INICIALIZA��O DAS VARI�VEIS DO PROBLEMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nelx=30; nely=15;
vol = 1; change = 1; ij = 0;
figx=10; figy=10; %%Quantidade de colunas e linhas a serem usadas no display da popula��o em cada gera��o
%% DADOS GEN�TICOS
rng('shuffle') 
itlim=100;
maxfitness(1:itlim)=0; %%Melhor fitness de cada itera��o
bestfitness(1:itlim)=0; %%Melhor fitness global
mutationrate=0.10; %%Porcentagem do genoma que pode ser alterado em cada itera��o
popsize=100; %%Tamanho da popula��o
chromosome=zeros(nelx*nely,popsize); %%Pr�-alocamento dos genes, a ser esvaziado posteriormente
deathrate=0.95; %%Quantos porcento da popula��o morre em cada gera��o - o complementar para 1 de quantos se mant�m
liven=round((1-deathrate)*popsize);
rep_probability=zeros(1,popsize); %%Porcentagem de um cromosomo ser escolhido para reprodu��o
passeddowngenes=0.6; %%Porcentagem de genes a serem passados do melhor pai para a crian�a
passeddowngenes=round(passeddowngenes*nely*nelx);
pairs=zeros(popsize,2); %%Array de rela��o de casais
%% DADOS DO MEIO MATERIAL
E0 = 1000; nu = 0.3; rho0 = 1000; %Dados de entrada
Lx = 2; Ly = 1; h = 0.01; %Dados do problema
xmin=1e-3; %penal=3; %%Par�metros do m�todo soft-kill
ctp=2; %%Determina as condi��es do contorno do meio
%% PREPARE FINITE ELEMENT ANALYSIS
a = 1/2*(Lx/nelx); b = 1/2*(Ly/nely);
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE = h/(1-nu^2)/24*([A11 A12; A12' A11]+nu*[B11 B12; B12' B11]);
ME = 1/9*a*b*h*toeplitz([4 0 2 0 1 0 2 0]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS
switch(ctp)
    case 1 % HALF-MBB BEAM
        F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
        fixeddofs = union(1:2:2*(nely+1),2*(nely+1)*(nelx+1));
    case 2 % CANTILEVER
        F = sparse(2*(nely+1)*(nelx+1)-nely,1,-1,2*(nely+1)*(nelx+1),1);
        fixeddofs = (1:2*(nely+1));
    case 3 % HALF-WHEEL
        F = sparse(2*(nely+1)*(nelx/2+1),1,-1,2*(nely+1)*(nelx+1),1);
        fixeddofs = union(2*nely+1:2*(nely+1),2*(nely+1)*(nelx+1));
    case 4 % OLHOFF BEAM
        F = sparse((1+nelx)*nely+nelx+2,1,-1,2*(nely+1)*(nelx+1),1);
        fixeddofs = union(nely+1:nely+2,(1+2*nelx)*nely+2*nelx+1:(1+2*nelx)*nely+2*nelx+2);
    case 5 % CLAMPED BEAM
        F = sparse((1+nelx)*nely+nelx+2,1,-1,2*(nely+1)*(nelx+1),1);
        fixeddofs = union(1:2*(nely+1), 1+2*(nely+1)*nelx:2*(nely+1)+2*(nely+1)*nelx);
    case 6 % A retangular Plate Steven 1996
        F = sparse((1+nelx)*nely+nelx+2,1,-1,2*(nely+1)*(nelx+1),1);
        fixeddofs =  union( 2*nely+1:2*nely+2 , (2*(nelx+1)*(nely+1)-2*nely-1):2*(nelx+1)*(nely+1)-2*nely );
end
U = zeros(2*(nely+1)*(nelx+1),1);
alldofs = (1:2*(nely+1)*(nelx+1));
freedofs = setdiff(alldofs,fixeddofs);

%% IN�CIO DO PROGRAMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INICIALIZAR A POPULA��O INICIAL
population=randi([0,1],nely,nelx,popsize);
%population(:,:,5)=ones(nely,nelx);
population(:,:,1)=zeros(nely,nelx);
population(floor(nely/2):ceil(nely/2),:,1)=ones(2,nelx);
while ij<itlim
    ij = ij + 1;
    %% AN�LISE DE FITNESS DE CADA GERA��O
    for individual=1:popsize
        chromosome(:,individual)=reshape(population(:,:,individual),nelx*nely,1);
        %%AN�LISE DE ELEMENTOS FINITOS 
        % VOLTADA PARA FREQU�NCIAS
        %     sK = reshape(KE(:)*((xmin-xmin^penal)/(1-xmin^penal)*(1 - x(:)'.^penal) + x(:)'.^penal)*E0,64*nelx*nely,1);
        %     sM = reshape(ME(:)*x(:)'*rho0,64*nelx*nely,1);
        %     K = sparse(iK,jK,sK); K = (K+K')/2;
        %     M = sparse(iK,jK,sM); M = (M+M')/2;
        %     if ctp == 5; M(nely+2+nelx*(nely+1),nely+2+nelx*(nely+1)) = mass; end
        % VOLTADA PARA COMPLIANCE/RIGIDEZ
        sK = reshape(KE(:)*max(xmin,chromosome(:,individual))'*E0,64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U(freedofs) = K(freedofs,freedofs)\F(freedofs);
        weight(individual)=rho0*mean(mean(population(:,:,individual)))*h;
        % AN�LISE DE COMPLIANCE DO INDIV�DUO
        ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
        distance(individual)=max(U(freedofs));
        c(ij,individual) = sum(sum((population(:,:,individual).*E0).*ce));
        prefit(ij,individual)=distance(individual)*weight(individual);
        % CHECAGEM DE CONEX�O ENTRE PONTOS
        connection=1;
        for i=1:nely
            for j=1:nelx
                if population(i,j)==1
                    if i~=1&&i~=nely
                        if j~=1&&j~=nelx
                            if population(i+1,j)==0&&population(i-1,j)==0&&population(i,j+1)==0&&population(i,j-1)==0
                                connection=1e-5;
                            end
                        end
                    end
                end
            end
        end
        prefit(ij,individual)=prefit(ij,individual).^-1*connection;
    end
    %%CALCULO DA FUN��O FITNESS
    % -Deve haver alguma maneira de tornar isso melhor (checando
    % conectividade, manufaturabilidade, etc), mas como eu sou pregui�oso
    % e isso � um projeto de fim de semana, � so stiffness mesmo,
    % meu br�di truta B^)
    
    %%Stifness=Compliance^-1 ======== Fitness
    [fitness,ordem]=sort(prefit(ij,:),'descend');
    maxfitness(ij)=fitness(1);
    if maxfitness(ij)>bestfitness(ij)
        bestfitness(ij)=maxfitness(ij);
    else
        bestfitness(ij)=bestfitness(ij-1);
    end
    %%DISPLAY DOS [FIGY*FIGX] INDIV�DUOS COM FITNESS
    fprintf('Generation:%3i Generation Fitness:%8.4f Best Fitness:%8.4f \n',(ij),maxfitness(ij),bestfitness(ij));
    for i=1:figy
        for j=1:figx
            subplot(figy,figx,figx*(i-1)+j);
            colormap(gray); imagesc(-population(:,:,ordem(figx*(i-1)+j))); axis equal; axis tight; axis off;
            %title(['Fitness: ' num2str(fitness(figx*(i-1)+j)) ]);
        end
    end
    population=population(:,:,ordem);
    chromosome=chromosome(:,ordem);
    pause(1e-6)
    %% SELE��O DOS INDIV�DUOS
    %%CALCULO DE PROBABILIDADES DE REPRODUCAO
    filter=sum(fitness); %%Fitness proportionate selection
%     filter=sum(1:popsize);
    for i=1:popsize
        rep_probability(i)=fitness(i)/filter;
        %rep_probability(i)=(popsize-i+1)/filter;
        cumprob(i)=sum(rep_probability(1:i));
    end
    
    %%SELE��O DOS PAIS
    for i=1:popsize
        lucky=rand;
        j=1;
         %%SELE��O DO PRIMEIRO PARENTE
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
        %%SELE��O DO SEGUNDO PARENTE
        % NORMAL, SEM PREVEN��O DE INCESTO
        lucky=rand;
        j=1;
        if lucky>=cumprob(j)
            while lucky>=cumprob(j)
                j=j+1;
                if lucky<cumprob(j) %%exato mesmo c�digo pq pregui�a
                    pairs(i,2)=j;
                    break
                end
            end
        else
            pairs(i,2)=j;
        end
        
        % PREVEN��O DE INCESTO CONSIGO MESMO (EXPERIMENTAL)
        % selmethod=randi([0,rep_probability(j)]);
        % tempprob=cumprob-rep_probability(j);
        
        % CHECAGEM DE DOMIN�NCIA E CORRE��O
        if fitness(pairs(i,2))>fitness(pairs(i,1))
            temp=pairs(i,2);
            pairs(i,2)=pairs(i,1);
            pairs(i,1)=temp;
        end
    end
    
    %% CROSSOVER DOS GENES
    oldpopulation=population;
    oldchromosome=chromosome;
    %%PAR RECESSIVO COMO BASE, SER� ACRESCENTADO OS GENES DOMINANTES
    for i=1:popsize
        crossover(:,i)=chromosome(:,pairs(i,2));
    end
    
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
    for i=1:popsize
        rng('shuffle')                                                         %Garantia de passeddowngenes%
        dominantgenes=datasample(1:nely*nelx,passeddowngenes,'replace',false); %dos genes dominantes serem passados
        for j=1:passeddowngenes 
            crossover(dominantgenes(j),i)=chromosome(dominantgenes(j),pairs(i,1));
        end
        population(:,:,i)=reshape(crossover(:,i),[nely,nelx]);
        chromosome(:,i)=crossover(:,i);
    end
    
    %% MUTATIONS
    %%DETEC��O SE EST� CONVERGINDO DEMAIS
%     if ij>10
%         if abs(sum(maxfitness(ij-9:ij-5))-sum(maxfitness(ij-4:ij)))/sum(maxfitness(ij-4:ij))<0.3
%             mutation=mutation+0.05;
%         else
%             mutation=mutation-0.05;
%         end
%     end
    %%MUTA��O BITFLIP SIMPLES
    for i=1:popsize
        mutation=rand;
        %%LEITURA DOS CROMOSOMOS
        for m=1:(nelx*nely)
            if mutation<mutationrate
                %%MUTA��O
                if chromosome(m,i)==1
                    chromosome(m,i)=0;
                else
                    chromosome(m,i)=1;
                end
            end
        end
        %%REFORMA��O
        population(:,:,i)=reshape(chromosome(:,i),[nely,nelx]);
    end

    %% OS PECADOS DOS PAIS
    for sin=1:liven
        population(:,:,(sin))=oldpopulation(:,:,sin);
    end
end

figure()
colormap(gray); imagesc(-population(:,:,1)); axis equal; axis tight; axis off;
title(['Max Fitness of the last generation: ' num2str(maxfitness(ij)) ' Pa'])
figure()
plot(1:ij,maxfitness(ij-1),'bo')
hold on
plot(1:ij,bestfitness(ij-1),'kd')
grid on
grid minor