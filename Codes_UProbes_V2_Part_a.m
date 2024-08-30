
% Date    |> 2024.05.13
% Project |> Non-complementary Strand Commutation
% Author  |> Guan Alex Wang
% 
% Version Info 
%       -> [1]. Function - 1: Typical Yield ~ dG_rxn for TER 
%       -> [2]. Function - 2: Calculatet the Yield ~ (at) for DNA Commutation
%       -> [3]. Function - 3: Calculatet the Yield ~ (at) for DNA Commutation and SNV
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

STD_dG_TER = linspace(-15, 15, 301)';
STD_K_TER = exp(-STD_dG_TER./(R*T_kel));

STD_simY_TER = zeros(size(STD_dG_TER) );

pa = struct('T_0', 12.5, 'CP_0', 1, 'P_0', 9); % using 10nM as reference

wb_TER = waitbar(0, 'Please wait...');
x = sym('x', [1, 1]);
assume(x, 'real')
x_lb = min(pa.T_0, 1);

for i = 1: numel(STD_dG_TER)
    Eqn_TER = ...
        x*(pa.P_0+x)/( (pa.T_0-x)*(pa.CP_0-x) ) == STD_K_TER(i);
    
    sol_TER = vpasolve(Eqn_TER, x, [0, x_lb] );
    STD_simY_TER(i) = sol_TER(sol_TER >= 0 & sol_TER < x_lb );
    
    waitbar(i/numel(STD_dG_TER), wb_TER, sprintf('Process: %d / %d',i,numel(STD_dG_TER)) )
end
close(wb_TER)

figure
plot(STD_dG_TER, STD_simY_TER, 'o', 'LineWidth', 1);
set(gca,'Xdir','reverse')


%%                              |>  Main - Typical NCSM  <|

% |> Association 1: C + P <-> CP, dG_A1, K = exp(-dG_A1/RT)*2e-8
% |> Association 2: C + T <-> CT, dG_A2, K = exp(-dG_A2/RT)*2e-8

mat_rxn = [-1 -1 0 1 0; 0 1 -1 -1 1; 1 0 1 0 -1];
rank(mat_rxn)

% ................................................................................

STD_dGa_CP = linspace(-15, -5, 101)';
STD_dGa_TC = linspace(-15, -5, 101)';
STD_Ka_CP = exp(-STD_dGa_CP./(R*T_kel)) * 2e-8;
STD_Ka_TC = exp(-STD_dGa_TC./(R*T_kel)) * 2e-8;

STD_simY_TC = zeros(length(STD_dGa_CP), length(STD_dGa_TC)); 
STD_simY_C = zeros(length(STD_dGa_CP), length(STD_dGa_TC)); 

y = sym('y', [2, 1]);
assume(y, 'real')
y_lb = min(pa.T_0, 1);

% |> setting: Y_TC = y1; Y_C = y2
for j = 1: length(STD_dGa_TC)
    for i = 1: length(STD_dGa_CP)
        Eqn_Asso_CP = (1-y(1)-y(2))/( y(2)*(pa.P_0+y(1)+y(2)) ) == STD_Ka_CP(i);
        Eqn_Asso_TC = y(1)/( y(2)*(pa.T_0-y(1)) ) == STD_Ka_TC(j);
        Eqn_TER = y(1)*(pa.P_0+y(1)+y(2))/((pa.T_0-y(1))*(1-y(1)-y(2))) == STD_Ka_TC(j)/STD_Ka_CP(i);

        sol_Asso = vpasolve([Eqn_Asso_CP, Eqn_Asso_TC], y, [0, y_lb; 0, y_lb] );
        temp_sol_y1 = sol_Asso.y1;
        temp_sol_y2 = sol_Asso.y2;

        STD_simY_TC(j, i) = temp_sol_y1(temp_sol_y1 >= 0 & temp_sol_y1 < y_lb);
        STD_simY_C(j, i) = temp_sol_y2(temp_sol_y2 >= 0 & temp_sol_y2 < y_lb);

        fprintf('Row: %d / %d - Column: %d / %d  \n', j,length(STD_dGa_TC), i,length(STD_dGa_CP) )
    end; clc
end

% ............................................................................................
[axis_X1, axis_Y1] = meshgrid(STD_dGa_CP, STD_dGa_TC);

STD_Ymat_TC_temp = STD_simY_TC;
STD_Ymat_C_temp = STD_simY_C;
STD_Ymat_FL_temp = STD_simY_TC + STD_simY_C;

STD_Ymat_TC_temp(STD_Ymat_TC_temp <= 0.01 ) = 0.01;
STD_Ymat_C_temp(STD_Ymat_C_temp <= 0.01 ) = 0.01;
STD_Ymat_FL_temp(STD_Ymat_FL_temp <= 0.01 ) = 0.01;

figure
m_simY_matTC = pcolor(axis_X1, axis_Y1, STD_Ymat_TC_temp);
c_simY_matTC = colorbar;
c_simY_matTC.Label.String = 'Sim Yield';
set(m_simY_matTC, 'LineStyle', 'none');
% set(gca,'Xscale','log');
% set(gca,'Yscale','log');
clim([0 1]);
ylabel('Association dG_{TC} / (kcal/mol)','FontSize',10);
xlabel('Association dG_{CP} / (kcal/mol)','FontSize',10);
title('Yield of TC','FontSize',12);
colormap('jet');

figure
m_simY_matC = pcolor(axis_X1, axis_Y1, STD_Ymat_C_temp);
c_simY_matC = colorbar;
c_simY_matC.Label.String = 'Sim Yield';
set(m_simY_matC, 'LineStyle', 'none');
% set(gca,'Xscale','log');
% set(gca,'Yscale','log');
clim([0 1]);
ylabel('Association dG_{TC} / (kcal/mol)','FontSize',10);
xlabel('Association dG_{CP} / (kcal/mol)','FontSize',10);
title('Yield of C');
colormap('jet');


% overall Fluorescence
figure
m_simY_matFL = pcolor(axis_X1, axis_Y1, STD_Ymat_FL_temp);
c_simY_matFL = colorbar;
c_simY_matFL.Label.String = 'Sim Yield';
set(m_simY_matFL, 'LineStyle', 'none');
% set(gca,'Xscale','log');
% set(gca,'Yscale','log');
clim([0 1]);
ylabel('Association dG_{TC} / (kcal/mol)','FontSize',10);
xlabel('Association dG_{CP} / (kcal/mol)','FontSize',10);
title('Yield of FL');
colormap('jet');


%% Sim Data Analysis

STD_simY_bkg = zeros(size(STD_dGa_CP))';

for i = 1: length(STD_dGa_CP)
    Eqn_Asso_CP = (1-x)/( x*(pa.P_0+x) ) == STD_Ka_CP(i);
    sol_bkg = vpasolve(Eqn_Asso_CP, x, [0, 1] );
    
    STD_simY_bkg(i) = sol_bkg;
    fprintf('Row: %d \n', i, length(STD_dGa_CP) );
end; clc
STD_Ymat_bkg_temp = STD_simY_bkg;
STD_Ymat_bkg_temp(STD_Ymat_bkg_temp <= 0.01 ) = 0.01;


% Signal / Noise ratio
figure
m_mat_SNr = pcolor(axis_X1, axis_Y1, STD_Ymat_FL_temp./STD_Ymat_bkg_temp);
c_mat_SNr = colorbar;
c_mat_SNr.Label.String = 'Sim S/N Ratio';
set(m_mat_SNr, 'LineStyle', 'none');
% set(gca,'Xscale','log');
% set(gca,'Yscale','log');
clim([1 10]);
ylabel('Association dG_{TC} / (kcal/mol)','FontSize',10);
xlabel('Association dG_{CP} / (kcal/mol)','FontSize',10);
title('Yield of S/N Ratio');
colormap('jet');


%% SNV Discrimination

ddG_snv = [1, 2, 3, 4, 5, 6, 7];

STD_dGa_TC_snv = STD_dGa_TC*ones(size(ddG_snv)) + ones(size(STD_dGa_TC))*ddG_snv;  % a Matrix
STD_Ka_TC_snv = exp(-STD_dGa_TC_snv./(R*T_kel)) * 2e-8;

STD_simY_TC_snv = zeros(length(STD_dGa_CP), length(STD_dGa_TC), length(ddG_snv));  % a 3D Matrix
STD_simY_C_snv = zeros(length(STD_dGa_CP), length(STD_dGa_TC), length(ddG_snv));  % a 3D Matrix

% |> setting: Y_TC = y1; Y_C = y2
for l = 1: length(ddG_snv)
    for j = 1: length(STD_dGa_TC)
        for i = 1: length(STD_dGa_CP)
            Eqn_Asso_CP = (1-y(1)-y(2))/( y(2)*(pa.P_0+y(1)+y(2)) ) == STD_Ka_CP(i);
            Eqn_Asso_TC_snv = y(1)/( y(2)*(pa.T_0-y(1)) ) == STD_Ka_TC_snv(j, l);
%             Eqn_TER = y(1)*(pa.P_0+y(1)+y(2))/((pa.T_0-y(1))*(1-y(1)-y(2))) == STD_Ka_TC(j)/STD_Ka_CP(i);

            sol_Asso_snv = vpasolve([Eqn_Asso_CP, Eqn_Asso_TC_snv], y, [0, y_lb; 0, y_lb] );
            temp_sol_y1 = sol_Asso_snv.y1;
            temp_sol_y2 = sol_Asso_snv.y2;

            STD_simY_TC_snv(j, i, l) = temp_sol_y1(temp_sol_y1 >= 0 & temp_sol_y1 < y_lb);
            STD_simY_C_snv(j, i, l) = temp_sol_y2(temp_sol_y2 >= 0 & temp_sol_y2 < y_lb);

            fprintf('SNV ID# : %d  |  Row: %d / %d - Column: %d / %d  \n', l, ...
                j,length(STD_dGa_TC), i,length(STD_dGa_CP) )
        end; clc
    end
end


%% Sim Data Analysis of SNV

STD_Ymat_TC_snv_temp = STD_simY_TC_snv;
STD_Ymat_C_snv_temp = STD_simY_C_snv;
STD_Ymat_FL_snv_temp = STD_simY_TC_snv + STD_simY_C_snv;

STD_Ymat_TC_snv_temp(STD_Ymat_TC_snv_temp <= 0.01 ) = 0.01;
STD_Ymat_C_snv_temp(STD_Ymat_C_snv_temp <= 0.01 ) = 0.01;
STD_Ymat_FL_snv_temp(STD_Ymat_FL_snv_temp <= 0.01 ) = 0.01;

% ............................................................................................
for l = 7: length(ddG_snv)
    fprintf('the ddG of SNV is:  %d  kcal/mol ; \n', ddG_snv(l) )

    figure
    m_simY_matTC_snv = pcolor(axis_X1, axis_Y1, STD_Ymat_TC_snv_temp(:,:,l));
    c_simY_matTC_snv = colorbar;
    c_simY_matTC_snv.Label.String = 'Sim Yield';
    set(m_simY_matTC_snv, 'LineStyle', 'none');
    % set(gca,'Xscale','log');
    % set(gca,'Yscale','log');
    clim([0 1]);
    ylabel('Association dG_{TC} / (kcal/mol)','FontSize',10);
    xlabel('Association dG_{CP} / (kcal/mol)','FontSize',10);
    title('Yield of TC of SNV','FontSize',12);
    colormap('jet');


    figure
    m_simY_matC_snv = pcolor(axis_X1, axis_Y1, STD_Ymat_C_snv_temp(:,:,l));
    c_simY_matC_snv = colorbar;
    c_simY_matC_snv.Label.String = 'Sim Yield';
    set(m_simY_matC_snv, 'LineStyle', 'none');
    % set(gca,'Xscale','log');
    % set(gca,'Yscale','log');
    clim([0 1]);
    ylabel('Association dG_{TC} / (kcal/mol)','FontSize',10);
    xlabel('Association dG_{CP} / (kcal/mol)','FontSize',10);
    title('Yield of C of SNV');
    colormap('jet');


    % overall Fluorescence
    figure
    m_simY_matFL_snv = pcolor(axis_X1, axis_Y1, STD_Ymat_FL_snv_temp(:,:,l));
    c_simY_matFL_snv = colorbar;
    c_simY_matFL_snv.Label.String = 'Sim Yield';
    set(m_simY_matFL_snv, 'LineStyle', 'none');
    % set(gca,'Xscale','log');
    % set(gca,'Yscale','log');
    clim([0 1]);
    ylabel('Association dG_{TC} / (kcal/mol)','FontSize',10);
    xlabel('Association dG_{CP} / (kcal/mol)','FontSize',10);
    title('Yield of FL of SNV');
    colormap('jet');


    % calculating DF values
    figure
    m_simY_matDF_snv = pcolor(axis_X1, axis_Y1, STD_Ymat_FL_temp./STD_Ymat_FL_snv_temp(:,:,l));
    c_simY_matDF_snv = colorbar;
    c_simY_matDF_snv.Label.String = 'Sim Yield';
    set(m_simY_matDF_snv, 'LineStyle', 'none');
    % set(gca,'Xscale','log');
    % set(gca,'Yscale','log');
    clim([1 20]);
    ylabel('Association dG_{TC} / (kcal/mol)','FontSize',10);
    xlabel('Association dG_{CP} / (kcal/mol)','FontSize',10);
    title('DF of SNV');
    colormap('jet');

end










