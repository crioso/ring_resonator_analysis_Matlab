%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file returns the ER, Q, FSR, ng, loss for under, over and critically
%%% coupled regimes, and the values for each pick Lorentz fitting. 
%%%
%%% How to use:
%%% 1. Check the global parameters section and enter, among others, the ring
%%%    radious, . 
%%% 2. Run the code and follow the instructions of the pop-up windows
%%%  - Fist pop-up window is to define the minimum and maximum wavelength
%%%    (i.e. to enter the wavelength range of interest). It plots the full
%%%    spectrum of the first device for easy visualization. This window also
%%%    allows to choose whether or not to do a Gaussian fit (in the case of
%%%    using grating couplers)
%%%  - Second pop-up window is to enter the approximate separation of
%%%    consecutive peask (i.e. ~FSR) to make sure only the real peaks are
%%%    caught as local minima of interest (sometime noise can be picked up
%%%    and their Q is huge...)
%%% 3. Go to the results folder and check the .xls file with all the
%%%    results and also the .jpg for the spectrum and each peak with the
%%%    fitting parameters. 
%%%
%%% All '*.dat' sheets in the current folder are imported. 
%%% The code is designed for files with wavelength in column 1 and
%%% tranmission data in column 2. It won't work otherwise. 
%%%
%%% Important! enter the wavelength step in global parameters. Usually 1pm.
%%%
%%% By: Carlos Ríos / MIT - April 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc


%% global parameter
linesheader = 1; %number of lineheaders to avoid from the data file
addpath('./lib'); %path for the library with the functions used in this code

radius = 120.0e+3;   %ring radius in nm
length = 2.0*pi*radius; %ring length in nm

wavelength_step=0.001; % wavelength step in the data collection process, given in nm

%% Use or creates the folder 'results' to store the data
folder_results_lorentz = 'results'; 
if isdir(folder_results_lorentz) == 0
    mkdir(folder_results_lorentz);
end
format long
%% Import data parameters
% order of the files during import does not matter since the results are
% stored with the original file name. Every file in the folder is imported!
testfiledir = '.\';
matfiles = dir(fullfile(testfiledir, '*.csv'));
%%%%%%%%%%% First guesses for the results -> used for the fits

B0_lorentz = 1;  % [nm] Used in fit_lorentz to determine the start parameters for the fit. This parameter is related to the FWHM and it's fine as long as it larger than the real FWHM so change if 1nm is too small. However, it's good to inspect the lorentz fittings, if the fitting is poor, try changing this parameter. 

%% Parameters for figures
scrsz = get(0,'ScreenSize');

label_x_axis = 'wavelength \lambda [nm]'; %label X for data subplots
label_y_axis = 'Power P [mW]'; %label Y for data subplots

fig1=figure('Name','Data_Gauss','Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/3 scrsz(4)/3]); %position and size of the subplots
fig2=figure('Name','Lorentz','Position',[scrsz(3)/2 scrsz(4)/4 scrsz(3)/3 scrsz(4)/3]);

%% Choose the wavelenght range
% Plots the first data file for a decision on wavelength
data_temp=csvread(matfiles(1).name,linesheader);  
    
figure;
plot(data_temp(:,1),10.^(data_temp(:,2)/10))
title('First device spectrum just to choose wavelength range')
xlabel(label_x_axis);
ylabel(label_y_axis);

prompt={'Enter minimum wavelength:','Enter maximum wavelength:','Gaussian Fit? (yes/no)'};
name='Wavelength range in nm';
numlines=1;
defaultanswer={'1540','1580','yes'};

answer=inputdlg(prompt,name,numlines,defaultanswer);
lambda_min=str2double(answer{1});
lambda_max=str2double(answer{2});
gaussian_cond=answer{3};
 
delete data_temp prompt name numlines defaultanswer answer
%% Import files and processing
[numberoffiles,aa]=size(matfiles)
for i=1:numberoffiles
        
    %Use/create folder to save the jpg for each peaks per device
    folder_temp = strcat(folder_results_lorentz,'/peak_fit_',matfiles(i).name);
    if isdir(folder_temp) == 0   
           mkdir(folder_temp);
    end
    
    %imports data
    data_temp=csvread(matfiles(i).name,linesheader);  
%     data_temp =flip(data_temp);   % In case the order has to be inverted
    data_temp_x0 = data_temp(:,1); % takes first column -> wavelength
    row_min = find(data_temp_x0>lambda_min-wavelength_step & data_temp_x0<lambda_min+wavelength_step/2); % Finds the place in the vector of the min wavelength
    row_max = find(data_temp_x0>lambda_max-wavelength_step & data_temp_x0<lambda_max+wavelength_step/2); % Finds the place in the vector of the max wavelength
    
    data_temp_x = data_temp(row_min:row_max,1); %Data for each device to be analyzed (wavelength)
    data_temp_y = 10.^(data_temp(row_min:row_max,2)/10); %Data for each device to be analyzed (Power in mW)
    data_temp_y_dB = data_temp(row_min:row_max,2); %Data for each device to be analyzed (Power in dB)
    numberdatapoints = size(data_temp_x,1);
    
    delete data_temp_x0 row_min row_max data_temp
    
    % Plots the transmission within the wavelength range of interest
    figure(fig1);
    plot(data_temp_x,data_temp_y,'r');
    xlabel(label_x_axis);
    ylabel(label_y_axis);
    hold('on');

    % Fit the spectrum with a Gaussian function
    if strcmp(gaussian_cond,'yes')
        curvefit = fit(data_temp_x,data_temp_y,'gauss1'); % see cflibhelp for description of gauss1
        x0 = curvefit.b1; % use gauss(x) = A*exp((x-x0)^2/sigma^2)
        sigma = curvefit.c1;
        A = curvefit.a1;
        Amax = max(data_temp_y);

        g=@(x) Amax*exp(-(x-x0)^2/sigma^2);

        %%%%%%% Plot the fit

        fplot(g,[data_temp_x(1),data_temp_x(numberdatapoints)],'b');
    %     fplot(g,data_temp);
        hold('off');
        spectrum_pic = strcat(folder_temp,'/spectrum',num2str(i),'.jpg');
        print('-djpeg',char(spectrum_pic));

        %%%%%%% Normalize the curve by the gaussfit
        data_mod = zeros(numberdatapoints,2);  % normalized curve

        for k=1:numberdatapoints
            data_mod(k,1) = data_temp_x(k);
            if g(data_temp_x(k))> (A*10^-3)                     % cut-off because the error increases far away from the center
                data_mod(k,2) = data_temp_y(k)/g(data_temp_x(k));
            else
                data_mod(k,2) = 0;
            end
        end

    else %%%% Skip the normalization
         data_mod(:,1) = data_temp_x;
         data_mod(:,2) = data_temp_y;
    end
    
    %%% plot the data
    figure(fig2);
    plot(data_mod(:,1),data_mod(:,2),'b');
    xlabel(label_x_axis);
    ylabel(label_y_axis);
  
    data_mod_x = data_mod(:,1);
    data_mod_y = data_mod(:,2);

    %%%%%% Plot the data with selected peaks and insert the value for Delta
    
    pass = 'y';
    while pass == 'y'
        clear min_x min_y n_min FSR_vec;
        plot(data_mod_x,data_mod_y,'r');
        xlabel('wavelength \lambda [nm]');
        ylabel('power P [mW]');
         hold('on');
        Delta = str2double(inputdlg('Enter an approx. peak separation in nm')); % This is important to avoid measuring noise as resonance peaks. 
%         Delta=0.6;
    
        %%%%%% Detect the minima and calculate Delta_lambda
    
%         [Max,Min] = peakdet(data_mod_y,Delta,data_mod_x);
        Min = islocalmin(data_mod_y,'MinSeparation',Delta*1000/2);
        if isempty(Min)
            waitfor(errordlg('No minima could be detected. Choose a smaller value!','modal'));
        else
            min_x=data_mod_x(Min); % x-values of the different minima
            min_y=data_mod_y(Min); % y-values of the different minima
            n_min = size(min_x,1);  % number of minima
            
            %%%%%% Calculate free spectra range: 
            %%%%%% resonance peak for: m*lambda_m = n*L     therefore: lambda_m > lambda_m+1
            %%%%%% FSR = lambda_m+1 - lambda_m+1 = (lambda_m)^2 / (n*L + lambda_m)
            FSR_vec = zeros(n_min-1,1);
            for j=1:(n_min-1)
                FSR_vec(j) = min_x(j+1)-min_x(j);
            end
            
            %%%%%% Calculate effectiv group index: n_group = ((lambda_m^2 / FSR) - lambda_m)/L 
            n_effgroup_vec = zeros(n_min-1,1);
            for j=1:(n_min-1)
                n_effgroup_vec(j) = ( ((min_x(j+1)^2)/FSR_vec(j)) - min_x(j+1) )/length;
%                 n_effgroup_vec(j) = min_x(j+1)^2/FSR_vec(j)/length;
            end           

            temp_FSR_mean=0;
            for j=1:(n_min-1)
                temp_FSR_mean =temp_FSR_mean + FSR_vec(j);
            end
            FSR_mean = temp_FSR_mean/(n_min-1);

            %%%%%% Plot the detected minima to check their correctness

            for j=1:(n_min)
                plot(min_x(j),min_y(j),'x','LineWidth',2);
            end
            hold('off');
            button = questdlg('Do you want to change the current value (see wrong peaks)?', '','yes','no (go to the next file)', 'yes');
            switch button
                case 'yes'
                     pass = 'y';
                case 'no (go to the next file)'
                     pass = 'n';
            end
        end 
    end
   
    n_FSR = find(data_mod_x>data_mod_x(1)+FSR_mean-wavelength_step & data_mod_x<data_mod_x(1)+FSR_mean+wavelength_step/2);

    %Definitions
    FWHM = zeros(n_min,1);
    ER= zeros(n_min,1);
    Re = zeros(n_min,1);
    Q = zeros(n_min,1);
    
    % Takes each of the peaks selected in the previous step and performs a
    % Lorentz Fit
    for j=1:n_min
        
        i_min = find(data_mod_x>min_x(j)-wavelength_step & data_mod_x<min_x(j)+wavelength_step/2);

        %%%%%% Check if there are enough data points to the left and to the right
        
        if i_min < n_FSR/2
            continue;
        end
        if (i_min + n_FSR/2) > size(data_mod_x,1)
            continue;
        end
         
        x_mod = data_mod_x(i_min-round((n_FSR)/2):i_min+round((n_FSR/2)));
        y_mod = data_mod_y(i_min-round((n_FSR)/2):i_min+round((n_FSR/2)));
        
        [A_lorentz,B_lorentz,C_lorentz,x0_lorentz,SSE_lorentz] = fit_lorentz([x_mod,y_mod],min_x(j),min_y(j),B0_lorentz,min_x(j)-FSR_mean/3,min_x(j)+FSR_mean/3,min_x(j)-FSR_mean/3,min_x(j)+FSR_mean/3);
        FWHM(j) = 2*B_lorentz;
        Re(j) = 10*log10(C_lorentz/min_y(j));
        Q(j) = x0_lorentz/FWHM(j);
        
        % Plots and saves each resonance peak with the Lorentz Fit
        xlabel(label_x_axis);
        ylabel('transmission[normalized]');
        title(strcat('A = ',num2str(A_lorentz),' B = ',num2str(B_lorentz), ' C = ',num2str(C_lorentz),' x0 = ',num2str(x0_lorentz),' SSE = ',num2str(SSE_lorentz)));
        output_pic = strcat(folder_temp,'/lorentz_',num2str(j),'.jpg');
        print('-djpeg',char(output_pic));
    end
    
    %Definitions
    alpha_um = zeros(n_min-1,1);
    alpha_dBcm = zeros(n_min-1,1);
    alpha_under_dBcm=zeros(n_min-1,1);
    alpha_over_dBcm=zeros(n_min-1,1);
    k_under=zeros(n_min-1,1);
    k_over=zeros(n_min-1,1);
    
    for j=1:(n_min-1)
        if Q(j+1) == 0
            continue;
        end
        % Over and under-coupled losses
        [alpha_under_dBcm(j),k_under(j)] = Resonator_Loss(Q(j), Re(j), FSR_vec(j), length/1000, min_x(j)/1000, 0)
        [alpha_over_dBcm(j),k_over(j)] = Resonator_Loss(Q(j), Re(j), FSR_vec(j), length/1000, min_x(j)/1000, 1)
        % Calculate propagation loss coefficient (includes propagation loss and bend loss and coupling losses in case of loaded Q):
        % alpha = (2*pi*n_group)/(lambda*Q)
        alpha_um(j) = 2.0*pi*n_effgroup_vec(j)/(min_x(j+1)*1e-3*Q(j+1));
        alpha_dBcm(j) = 10/log(10)*1e4*alpha_um(j);
    end 
    
    % Saves all the results in a .csv file
    file_temp = char(strcat(folder_results_lorentz,'/results_',matfiles(i).name,'.csv'));
    [num_peaks,b]=size(Q)  
    varNames = {'Lambda','Q','ER', 'FSR','n_g','loaded loss (dBcm)', 'loss_under (dBcm)','k_under','loss_over (dBcm)','k_over'};
    T=table(min_x(1:num_peaks-1),Q(1:num_peaks-1),Re(1:num_peaks-1),FSR_vec,n_effgroup_vec(1:num_peaks-1),alpha_dBcm(1:num_peaks-1),alpha_under_dBcm,k_under,alpha_over_dBcm,k_over,'VariableNames',varNames);
    writetable(T,file_temp);
end


    