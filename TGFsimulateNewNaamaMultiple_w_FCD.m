function TGFsimulateNewNaamaMultiple_w_FCD(X,protein)
 
close all

%initialize structures and figures
BB = struct();
CC = struct();
DD = struct();
EE = struct();
CYTO = struct();
CYTOcomp = struct();

%load conditions
c = feval('TGFconditions_w_FCD');     
save conditions_beta.dat c -ascii;
%load parameters
p = feval('TGFparameters_w_FCD');
save parameters_beta.dat p -ascii;
pp=p;

%Dynamics values
tn    = c(1);     % Time span for integration,seconds
Tgfoff = c(2);
Tgfoff = 0.01;
Tgfbasal = c(3);
Tgfon = c(4);
tspan = [0:100:tn];
% tspan = [0;tn];
Tgfz = 's';

number_of_doses=2;
    if number_of_doses == 1
        TgfF = 1;
    else
        TgfF = log10(logspace(0.02,1,number_of_doses));
    end




iterations = 100;
ytc =1;
ztc = 1;


protlevelset = 1;
s3leveli=1;
    
    s3level = protlevelset(s3leveli);
        s4leveli=length(protlevelset);
    s4level = protlevelset(s4leveli);
    
    

s3varset = [0.1 0.2 0.4];
pvarset = [0.1 0.2 0.4 0.8];
% doseset = [0.2 0.5 1];
doseset = 0.2;
for dosei = 1:length(doseset)
%     Tgfon = TgfF(FFF+1);
    Tgfon = doseset(dosei);
    
        for s3vari = 1:length(s3varset);
            for pvari = 1:length(pvarset)
            s3var = s3varset(s3vari);   
            setvar = pvarset(pvari);


            % Tgfon = lognrnd(0,0.5,1,1).*1;

            %==========================================================================
            % Computing the unperturbed and perturbed solutions (dimensional solution)
            % Increasing one parameter at a time by X-fold


            S2nuc = cell(1,iterations); 
            S2cyto = cell(1,iterations); 
            S24cyto = cell(1,iterations); 
            S24nuc = cell(1,iterations); 
            Z = cell(1,iterations); 
            perturbationStrength = cell(1,iterations); 
            totalSmadz = cell(1,iterations);
            TimeVec = cell(1,iterations); 
            Params = cell(1,iterations);
            Dosez = cell(1,iterations); 
            Variationz = cell(1,iterations);


                parfor i=1:iterations;        % 11 parameters

                % Tgfon = lognrnd(0,1,1,1).*1;
                % Tgfon = unifrnd(0,1,1,1);
                disp(i) 

                %make parameter perturbation
                p = pp;
                %adjust the total protein level of Smad3 and Smad4
                p(9) = p(9).*s3level; %p(9) is smad3 total level. 
                p(10) = p(10).*s4level;

                variation = lognrnd(0,setvar,length(p),1);

                %     p(1:11)= pp(1:11); %only vary the fold-change parameters

                %vary all parameters and Smad3 by a specific amount
                %     variation = lognrnd(0,setvar,length(p),1); 
                variation(9) = lognrnd(0,s3var,1,1);
                p = p.*variation';
                % 

                p(15)=p(15).*ytc;
                p(17)=p(17).*ytc;
                p(18)=p(18).*ztc;
                p(16)=p(16).*ztc;




                %============
                %Time course for basal state
                %============
                y0 = TGFconcentrations_w_FCD(p);

                y0(22) = Tgfoff;
                %Computing initial guess for dy, using decic
                fixed_y0 = ones(size(y0));
                fixed_dy0 = zeros(size(y0));
                dy0 = zeros(size(y0));
                %                     [y0mod,dy0mod] = decic('TGFequations_w_FCD_nucratio',0,y0,fixed_y0,dy0,fixed_dy0);
                [y0mod,dy0mod] = decic(@(t,y,dy) TGFequations_w_FCD_nucratio(t,y,dy,p),0,y0,fixed_y0,dy0,fixed_dy0);
                % Solving the ODEs
                %                     [TT,YY] = ode15i('TGFequations_w_FCD_nucratio',tspan,y0mod,dy0mod);
                [TT,YY] = ode15i(@(t,y,dy) TGFequations_w_FCD_nucratio(t,y,dy,p),tspan,y0mod,dy0mod);

                %============
                %Time course for stimulated state
                %============    
                y0 = YY(end,:);    
                y0(22) = Tgfon;
                %  Computing initial guess for dy, using decic
                fixed_y0 = ones(size(y0));
                fixed_dy0 = zeros(size(y0));
                dy0 = zeros(size(y0));
                %                     [y0mod,dy0mod] = decic('TGFequations_w_FCD_nucratio',0,y0,fixed_y0,dy0,fixed_dy0);
                [y0mod,dy0mod] = decic(@(t,y,dy) TGFequations_w_FCD_nucratio(t,y,dy,p),0,y0,fixed_y0,dy0,fixed_dy0);
                % Solving the ODEs
                %                     [TTT,YYY] = ode15i('TGFequations_w_FCD_nucratio',tspan,y0mod,dy0mod);
                [TTT,YYY] = ode15i(@(t,y,dy) TGFequations_w_FCD_nucratio(t,y,dy,p),tspan,y0mod,dy0mod);

                %============
                %Concatentate the time courses for both states
                %============    
                T = vertcat(TT,TTT+TT(end));
                Y = vertcat(YY,YYY);



                %                 t=0;
                %                 y=0;
                %         
                %                 %============
                %                 %Time course for basal state
                %                 %============
                %                     yS0 = TGFconcentrations_w_FCD(p);
                %                     yS0(22) = Tgfoff;
                %                     % Solving the ODEs
                %              [TTs,YYs] = ode15s(@(t,y) TGFequations_w_FCD_nucratioz(t,y,p),tspan,yS0);
                % 
                %                 %============
                %                 %Time course for stimulated state
                %                 %============    
                %                     yS0 = YYs(end,:);    
                %                     yS0(22) = Tgfon;
                %                  % Solving the ODEs
                %                     [TTTs,YYYs] = ode15s(@(t,y) TGFequations_w_FCD_nucratioz(t,y,p),tspan,yS0);
                % 
                %                 %============
                %                 %Concatentate the time courses for both states
                %                 %============    
                %                 Ts = vertcat(TTs,TTTs+TTs(end));
                %                 Ys = vertcat(YYs,YYYs);               
                %                 
                %            








                totalSmad = p(9);
                % BB(i).Color = COLORS{perturbedParameter};


                %% species  
                details = chooseSpecies(T,Y);
                totalspecies = details.S2total;
                basal = find(T<max(tspan),1,'last');
                species = details.S2nuc;

                S2nuc{i} = details.S2nuc;
                S2cyto{i} = details.S2cyto;
                S24cyto{i} = details.S24cyto;
                S24nuc{i} = details.S24nuc;
                Z{i} = details.Z;
                perturbationStrength{i} = 10.^sum(abs(log10(variation)));
                totalSmadz{i} = p(9);
                TimeVec{i} = T;
                Params{i} = p;
                Dosez{i} = Tgfon;
                Variationz{i} = variation;





                % descriptors
                %                 CC = datastructmaker(protein,T,Y,basal,totalspecies,i,CC);
                %                 DD = datastructmaker('S24nuc',T,Y,basal,totalspecies,i,DD);
                %                 EE = datastructmaker('Z',T,Y,basal,totalspecies,i,EE);
                %                 CYTO = datastructmaker('S2cyto',T,Y,basal,totalspecies,i,CYTO);
                %                 CYTOcomp = datastructmaker('S24cyto',T,Y,basal,totalspecies,i,CYTOcomp);
                % FF(i,:) = Y(:,25);

                stophere=1;
                end
            % figures    
            % figure(secondfigure)
            % scatter(vertcat(BB.PerturbationStrength),vertcat(CC.peak)./nanmedian(vertcat(CC.peak)));hold on
            % scatter(vertcat(BB.PerturbationStrength),vertcat(CC.foldchange)./nanmedian(vertcat(CC.foldchange)));hold on
            % scatter(vertcat(BB.PerturbationStrength),vertcat(CC.percen)./nanmedian(vertcat(CC.percen)));hold on
            % scatter3(vertcat(BB.Dose),vertcat(DD.percen),vertcat(CC.percen));hold on



        figure(991239)
        % plotSensitivityAnalysis(CC,p)
        % plotSensitivityAnalysis(BB,p)
        stophere=1;
        % save('/Users/frick/Documents/Goentoro_Lab/DATA/Modeling/2015_08_24 Supplement Modeling/Naama p vary 01 while S3 vary 04/CONDITIONS/DataAfterTGFsimulateNEW.mat')
        % save('/Users/frick/Documents/Goentoro_Lab/Writing/Information Paper/Submission to PNAS/Rebuttal to PNAS/FIGURES/Reviewer1 point1/w FCD/DataAfterTGFsimulateFCD_NEW_noiseNONROBUST.mat')
        % save(strcat('/Users/frick/Documents/Goentoro_Lab/Writing/Information Paper/Submission to PNAS/Rebuttal to PNAS/FIGURES/Reviewer3 point2/data/fcdcorrelations.mat'));

        save(strcat('/Users/frick/Documents/Goentoro_Lab/Writing/Information Paper/Submission to PNAS/Rebuttal to PNAS/FIGURES/Reviewer1 point2/NucCytoCorr_sims/data/nuccytocorrelations-ytc',num2str(ytc),'-ztc',num2str(ztc),'-s3exp',num2str(s3leveli),'-s4exp',num2str(s4leveli),'-s3noise',num2str(s3vari),'-pnoise',num2str(pvari),'-dose',num2str(dosei),'.mat'));
            end
        end

end
% for i=1:length(TgfF)
% IMat = INFOyo{i};
% for j = 1:size(IMat,1)
% subplot(3,3,j);scatter(randi(1000,[1 1000]),IMat(j,:));hold on
% end
% end

end



function details = chooseSpecies(T,Y)
details.S2nuc = Y(:,11)+Y(:,13)+Y(:,16)+Y(:,18); %nuclear Smad2
details.S2cyto =  Y(:,1)+Y(:,3)+Y(:,6)+Y(:,8); %cytoplasmic Smad2
details.S2total = Y(:,11)+Y(:,13)+Y(:,16)+Y(:,18)+Y(:,1)+Y(:,3)+Y(:,6)+Y(:,8); %total smad2

details.S4nuc = Y(:,15)+Y(:,16); %nuclear S4
details.S4cyto = Y(:,5)+Y(:,6); %cytoplasmic S4
details.S4total = Y(:,15)+Y(:,16)+Y(:,5)+Y(:,6); %total S4

details.S24nuc = Y(:,16); %S24 nuclear
details.S24cyto = Y(:,6);
details.S24total = Y(:,16)+Y(:,6);
details.Z = Y(:,25); % Z from FCD circuit
  end
    

  
function CC = datastructmaker(protein,T,Y,basal,totalspecies,i,CC)
details = chooseSpecies(T,Y);
species = details.(protein);
rate = gradient(species);% rate
CC(i).maxrate = max(rate(basal+1:end));%max rate
CC(i).maxrelrate = max(rate(basal+1:end))./species(basal);% relative rate
CC(i).foldchange = species(length(species))./species(basal);% fold change
CC(i).basilico = species(basal);% basal
CC(i).peak = species(length(species));% peak
CC(i).percen = CC(i).foldchange-1;%percent
CC(i).NT = species(length(species))./totalspecies(length(totalspecies));% nuclear/total
CC(i).maxpeak = max(species(basal+1:end));
CC(i).maxFC = max(species(basal+1:end))./species(basal);
end

    
