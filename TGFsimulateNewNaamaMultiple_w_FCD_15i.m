function TGFsimulateNewNaamaMultiple_w_FCD_15i(X,protein)
%try to vary import and export rate differently.



%set the directory from which to load other files
    mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
        if isempty(b)
            [~,b] = regexp(mdir,'\');
        end
    parentdir = mdir(1:b(end-1));

%initialize structures and figures
    BB = struct();
    CC = struct();
    DD = struct();
    EE = struct();

%load conditions
    c = feval('TGFconditions_w_FCD');     

%load parameters
    p = feval('TGFparameters_w_FCD');
    pp=p;

%Dynamics values
    tn    = c(1);     % Time span for integration,seconds
    Tgfoff = c(2);
    Tgfbasal = c(3);
    Tgfon = c(4);
    % tspan = [0:100:tn];
    tspan = [0;tn];
    Tgfz = 's';



    

ztcset = 1;
s3varset = [0.1 0.4];
setvarset = [0 0.0125/2 0.0125 0.025 0.05 0.1 0.2 0.3 0.4];
% INEXvarset = [0.05 0.1 0.2 0.4];
INEXvarset = 1;
doseset = [0.2 0.5 1];

inexrun = 1:length(INEXvarset);
setvarun = 1:length(setvarset);
% setvarun = 2;
s3varun =  1; %1:length(s3varset);
doserun = [1 3]; %1:length(doseset)
iterations = 200;
Tgfoff = 0.01;

protlevelset = 1;
    for setvari =setvarun
        for s3vari = s3varun
            for INEXvari = inexrun;
                for ztcvar = 1:length(ztcset);
                    for s4leveli = length(protlevelset)
                        for s3leveli = length(protlevelset)
                            for dosei = doserun
                                Tgfon = doseset(dosei);
                                s3level = protlevelset(s3leveli);
                                s4level = protlevelset(s4leveli);
                                
                                setvar = setvarset(setvari);
                                s3var = s3varset(s3vari);
                                
                                INEXvar=INEXvarset(INEXvari);%import and export parameters
                                ztc=ztcset(ztcvar);

                                %initialize variables
                                
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

                                    
                                %initialize smad3 expression profile variables
        %                         mixtime=10; %50
                                    parfor i=1:iterations;        % 11 parameters
                                        p=pp;
                                        
                                        disp(i) 

                                        %make parameter perturbation
                                        %adjust the total protein level of Smad3 and Smad4
                                        p(9) = p(9).*s3level;
                                        p(10) = p(10).*s4level;

                                        %vary non-FCD parameters
                                        variation = lognrnd(0,setvar,11,1); %set variation to only affect pathway parameters and not FCD
%                                         variation = lognrnd(0,setvar,length(p),1);
                                        variation(9) = lognrnd(0,s3var,1,1);
                                        
                                        
%                                         variation(1:2) = lognrnd(0,INEXvar,2,1); %vary import and export rates 
                                         
                                        
                                        %p = p.*variation';%no noise in the FCD parameters
                                        p(1:11) = p(1:11).*variation(1:11)';%no noise in the FCD parameters 

                                        %alter time scales of y and z in FCD
                                        p(15)=p(15).*INEXvar;
                                        p(17)=p(17).*INEXvar;
                                        p(18)=p(18).*ztc;
                                        p(16)=p(16).*ztc;


                                        %==========================================================================
                                        % Computing the unperturbed and perturbed solutions (dimensional solution)
                                        %==========================================================================
                                        %============
                                        %Time course for basal state
                                        %============
                                        y0 = TGFconcentrations_w_FCD(p);
                                        y0(22) = Tgfoff;
                                        fixed_y0 = ones(size(y0));
                                        fixed_dy0 = zeros(size(y0));
                                        dy0 = zeros(size(y0));

                                        % Solving the ODEs
                                        % [TT,YY] = ode15s(@(t,y) TGFequations_w_FCD_nucratio_15s(t,y,p,tcnoise,tcnt),tspan,y0);
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
                                        [TTT,YYY] = ode15i(@(t,y,dy) TGFequations_w_FCD_nucratio(t,y,dy,p),tspan+tspan(end),y0mod,dy0mod);
                                        %[TTT,YYY] = ode15s(@(t,y) TGFequations_w_FCD_nucratio_15s(t,y,p,tcnoise,tcnt),tspan+tspan(end),y0);

                                        %============
                                        %Concatentate the time courses for both states
                                        %============    
                                        T = vertcat(TT(1:end-1),TTT);
                                        Y = vertcat(YY(1:end-1,:),YYY);

                        %                 figure, plot(T,sum(Y(:,[1 3 6 8 8]),2)+1/2.3*sum(Y(:,[11 13 16 18 18]),2))


                                        details = chooseSpecies(T,Y);
                                        species = details.S2nuc;
%                                         figure(91593)
%                                         subplot(2,2,1);plot((T-tn(end))./60,species);hold on
%                                         xlim([-120 600])
%                                         subplot(2,2,2);plot((T-tn(end))./60,Y(:,25));hold on
%                                             ylim([0 0.2])
%                                             xlim([-120 600])



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

                                        stophere=1;
                                    end




                            save(strcat(parentdir,'data/nuccytocorrelations-INEXvar',num2str(INEXvari),'-ztc',num2str(ztc),'-s3exp',num2str(s3leveli),'-s4exp',num2str(s4leveli),'-s3noise',num2str(s3vari),'-pnoise',num2str(setvari),'-dose',num2str(dosei),'.mat'));
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
            end
        end
    end
end




function details = chooseSpecies(~,Y)
% details.S2nuc = Y(:,11)+Y(:,13)+Y(:,16)+Y(:,18); %nuclear Smad2
details.S2nuc = sum(Y(:,[11 13 16 18 18]),2); %nuclear Smad2
details.S2cyto =  sum(Y(:,[1 3 6 8 8]),2); %cytoplasmic Smad2
details.S2total = sum(Y(:,[1 3 6 8 8]),2)+1/2.3*sum(Y(:,[11 13 16 18 18]),2); %total smad2

details.S4nuc = Y(:,15)+Y(:,16); %nuclear S4
details.S4cyto = Y(:,5)+Y(:,6); %cytoplasmic S4
details.S4total = Y(:,15)+Y(:,16)+Y(:,5)+Y(:,6); %total S4

details.S24nuc = Y(:,16); %S24 nuclear
details.S24cyto = Y(:,6);
details.S24total = Y(:,16)+Y(:,6);
details.Z = Y(:,25); % Z from FCD circuit
  end
  
    
  
function CC = datastructmaker(protein,T,Y,basal,totalspecies,i,CC,basalidx)
details = chooseSpecies(T,Y);
species = details.(protein);
rate = gradient(species);% rate
CC(i).maxrate = max(rate(basal+1:end));%max rate
CC(i).maxrelrate = max(rate(basal+1:end))./species(basal);% relative rate

% CC(i).foldchange = species(length(species))./species(basal);% fold change
CC(i).foldchange = nanmean(species(end-length(basalidx):end))./nanmean(species(basalidx));

% CC(i).basilico = species(basal);% basal
CC(i).basilico = nanmean(species(basalidx));

% CC(i).peak = species(length(species));% peak
CC(i).peak = nanmean(species(end-length(basalidx):end));

CC(i).percen = CC(i).foldchange-1;%percent

CC(i).NT = (CC(i).peak)./totalspecies(length(totalspecies));% nuclear/total

maximumidx = find(species == max(species(basal+1:end)),1,'last');

frontidx = maximumidx - round(length(basalidx)./4);
backidx  = maximumidx + round(length(basalidx)./4);

if backidx > length(species)
    backidx = length(species);
end
peakidx = frontidx:backidx;

CC(i).maxpeak = nanmean(species(peakidx));
% CC(i).maxpeak = max(species(basal+1:end));
CC(i).maxFC = CC(i).maxpeak./CC(i).basilico;
CC(i).peakidx = peakidx;
stophere=1;
end



function [tcnt,tcnoise] = tcnoisefun(mixtime,Mean2,Sig2,dt,tlength)
    tau = (60*60*mixtime)./(-log(0.5)); %50 hours is a mixing time of ~35 hours
    mu2 = log((Mean2^2)./sqrt(Sig2.^2+Mean2^2));
    sigma2 = sqrt(log(Sig2.^2/(Mean2^2)+1));
    Tnoise = 0:dt:tlength;
    N = length(Tnoise);
    Var2 = normrnd(0,1,[1,N]); %(lambda,[1,N]);
    Val = zeros(1,N);
    Val(1) = Var2(1);
    f = exp(-dt/tau);
    
    for ii = 2:N;
        Val(ii) = Val(ii-1)*f+sqrt(1-f^2)*Var2(ii);
    end
    
    Val = exp(mu2+sigma2*Val);
%     figure, plot(Tnoise,Val)
%     mean(Val)
%     std(Val)
    tcnoise = Val;
    tcnt = Tnoise;

end