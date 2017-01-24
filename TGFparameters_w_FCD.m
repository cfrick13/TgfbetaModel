function p = TGFparameters_w_FCD
%

%==========================================================================
% The best-fit parameters in Schmierer et al., 2008
%==========================================================================

ytimescale = 1/5000;
ztimescale = 1/10000;
%%%%%%%%%%%%%%%%%%%%%
%test set of parameters
%%%%%%%%%%%%%%%%%%%%%%
p(1) = 0.0056;       % /s, kex, page 1 of Supp 0.0056;  
p(2) = 0.0026;       % /s, kin, page 1 of Supp 0.0026
p(3) = 0.000404;     % /nM-s, kphosp, Fig S4kex
p(4) = (0.016/8.7);    % /nM-s, kon = koff / Kdiss, Kdiss is in Fig S4
p(5) = 0.016;        % /s, koff, page 2 of Supp
p(6) = 5.7;          % CIF, Fig S4
p(7) = 0.00657;      % /nM-s, kdephosp, Fig S4
%increasing dephos rate increases robustness!!!
%and decreasesing decreases robustness
p(8) = 1;            % nM, PPase, Table S1 
p(9)= 89.1;             % nM, S2total           
p(10)= 101.6;        % nM, S4total       
p(11) = 0.074;          %/nM-s, kTGFbeta

p(12) = 10;     % kxy
p(13) = 10;     % kxz
p(14) = 10;     % kyz
p(15) = 1;     % betaY
p(16) = 2;     % betaZ
p(17) = 1;     % gammmaY;
p(18) = 100;     % gammaZ;
% p(19) = 1; %differential binding
p(20) = 0.1; %leaky transcription of Y
p(21) = 0.1; %leaky transcription of Z

p(15) = p(15).*ytimescale;
p(17) = p(17).*ytimescale;

p(16) = p(16).*ztimescale;
p(18) = p(18).*ztimescale;
p(19) = 1; %initial Y conc
p(20) = 1; %initial Z conc




% % %absolute change robustness
% % p = zeros(11,1);
% % 
% p(1) = 0.056;  
% p(2) = 0.000026;
% p(3) = 0.0000404; 
% p(4) = 0.016/8.7;    % /nM-s, kon = koff / Kdiss, Kdiss is in Fig S4
% p(5) = 0.00016;
% p(6) = 5.7;          % CIF, Fig S4
% p(7) = 0.0000657;
% p(8) = 0.001;
% p(9)= 89.1;             % nM, S2total
% p(10)= 101.6;        % nM, S4total
% p(11) = 0.074;          %/nM-s, kTGFbeta

% 
% 
% % %%%%%%%% PERTURBED!!!%%%%%%%%%%%%%
% % ABSOLUTE CONCENTRATION ROBUSTNESS IS ACHIEVED AT HIGH DOSES OF TGFBETA
% % AND DECREASING THE OFF RATE OF SMAD COMPLEX BINDING
% %%%%%%%
% % AND INCREASING THE CONCENTRATION OF SMADS
% %%%%%%%
% 
% 
% % p = zeros(11,1);
% p(1) = 0.0056;       % /s, kex, page 1 of Supp 0.0056;  
% % p(1) = 0.056;       % /s, kex, page 1 of Supp 
% p(2) = 0.0026;       % /s, kin, page 1 of Supp 0.0026
% % p(2) = 0.0007;       % /s, kin, page 1 of Supp
% p(3) = 0.000404;     % /nM-s, kphosp, Fig S4kex p(3) = 0.000404; 
% p(4) = 0.016/8.7;    % /nM-s, kon = koff / Kdiss, Kdiss is in Fig S4
% p(5) = 0.016;        % /s, koff, page 2 of Supp
% % p(6) = 5.7;          % CIF, Fig S4
% p(6) = 1;          % CIF, Fig S4
% p(7) = 0.00657;      % /nM-s, kdephosp, Fig S4
% p(8) = 1;            % nM, PPase, Table S1 
% % p(9)= 89.1;             % nM, S2total           
% % p(10)= 101.6;        % nM, S4total       
% % p(11) = 0.074;          %/nM-s, kTGFbeta
% p(9)= 200.1;             % nM, S2total           
% p(10)= 200.6;        % nM, S4total       
% p(11) = 0.074;          %/nM-s, kTGFbeta
