
% Date    |> 2024.04.15
% Project |> Under Complementary Probe (U-Probe)
% Version |> 2.0
% Author  |> Guan Alex Wang
% 
% Info 
%   -> [1]. Function - 1: TER Yield ~ (dG of TER, ddG of SNV) 

% -----------------------------------------------------------------------------------------------
% 
% NUPACK setting: 1M Na+, 0M Mg++; for Toehold Exchange Reaction
% 
% Experimental Settings:
% CP annealing: 20 nM C + 200 nM P -|> 20 nM CP + 180 nM free P
% typical concentration of T in Systematic Examination: 250 nM


%%                                    |>  Input  <|
clear; clc
R = 1.987*10^-3; % kcal*K^-1*mol^-1
T_kel = 273.15 + 25; % temperatue / K

%%                              |>  Main - Typical TER  <|
% Experimental Conditions
%   -> using 20 nM of CP as Concentration Reference
%   -> t: 250/20 = 12.5; p: 180/20 = 9;

STD_dG_TER = linspace(-10, 10, 51)';
STD_ddG_SNV = linspace(0, 5, 51)';

STD_dG_Rxn = ones(length(STD_ddG_SNV), 1)*STD_dG_TER' + STD_ddG_SNV*ones(1, length(STD_dG_TER));
STD_K_TER = exp(-STD_dG_Rxn./(R*T_kel));
 
STD_simY_TER = zeros(length(STD_ddG_SNV), length(STD_dG_TER) );
[xx, yy] = meshgrid(STD_dG_TER, STD_ddG_SNV);

pa = struct('T_0', 12.5, 'CP_0', 1, 'P_0', 9); % using 20nM as reference

wb_TER = waitbar(0, 'Please wait...');
x = sym('x', [1, 1]); % shape of the Symbol x: 1 x 1 array
assume(x, 'real')
x_lb = min(pa.T_0, 1);

tic
for i = 1: numel(STD_simY_TER)
    Eqn_TER = ...
        x*(pa.P_0+x)/( (pa.T_0-x)*(pa.CP_0-x) ) == STD_K_TER(i);
    
    sol_TER = vpasolve(Eqn_TER, x, [0, x_lb] );
    STD_simY_TER(i) = sol_TER(sol_TER >= 0 & sol_TER < x_lb );
    
    waitbar(i/numel(STD_simY_TER), wb_TER, ...
        sprintf('Process: %d / %d',i,numel(STD_simY_TER)) )
end
close(wb_TER); toc
% ..................................................................................................

STD_simY_TER_f = STD_simY_TER;
STD_simY_TER_f(STD_simY_TER_f <= 0.01) = 0.01;
% Figure (1) Y ~ (dG, ddG) -> the Yield distribution of TER
figure
Figure_1 = pcolor(xx, yy, STD_simY_TER_f);
Fig_1_c = colorbar;
Fig_1_c.Label.String = 'Sim Yield';
set(Figure_1, 'LineStyle', 'none');
clim([0.01 1]);
colormap('jet');
xlabel('dG / kcal/mol')
ylabel('ddG / kcal/mol' )
set(gca, 'XDir', 'reverse');


STD_simDY_TER_f = STD_simY_TER_f(1, :) - STD_simY_TER_f;
% Figure (2) dY ~ (dG, ddG) -> the Yield Difference distribution of TER
figure
Figure_2 = pcolor(xx, yy, STD_simDY_TER_f);
Fig_2_c = colorbar;
Fig_2_c.Label.String = 'Sim Yield';
set(Figure_2, 'LineStyle', 'none');
clim([0 1]);
colormap('jet');
xlabel('dG / kcal/mol')
ylabel('ddG / kcal/mol' )
set(gca, 'XDir', 'reverse');


STD_simDYr_TER_f = STD_simDY_TER_f./STD_simY_TER_f(1, :);
% Figure (3) dY/Y ~ (dG, ddG) -> the Yield Difference Ratio distribution of TER
figure
Figure_3 = pcolor(xx, yy, STD_simDYr_TER_f);
Fig_3_c = colorbar;
Fig_3_c.Label.String = 'Sim Yield';
set(Figure_3, 'LineStyle', 'none');
clim([0 1]);
colormap('jet');
xlabel('dG / kcal/mol')
ylabel('ddG / kcal/mol' )
set(gca, 'XDir', 'reverse');




























