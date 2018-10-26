clear all
close all
clc

rng default

%% Input folder
% Set the directory of the simulation
directory       = 'E:\Westerschelde_f\Domein2\0024\';
run_timestep    = 25920; % Number of timesteps for one loop
%% Manual input for the shoal margin collapse
Ec=0.8; % Shape of the shoal margin collapse, given by the eccentricity
Ec2=0.25; % Shape of the shoal margin collapse deposit, given by the eccentricity
psi=-0.05;
D50=2*10^-4; %Grain size Western Scheldt
Fmud=1; %Mud fraction (low(1/3), avg(1), high(3))
FSavg=5.3; % average flow slides a year in the Western Scheldt
Lsm=300; % total shoal margin length Western Scheldt
GridSize=300; %distance in meters
GridSize2=0:20:2000; %distance between deposit location and shoal margin

dt = 0.25;  % Timestep in D3D-computation

% Model specific ID name
try
    % Determine the name of the simulation
    fid             = fopen(strcat(directory, 'source\runid1'),'r');
    a               = fread(fid);
    ID1             = 'WS2';
    fclose (fid);
catch
    disp('error runID')
end

%% Run Delft 3D

for i = 1:100 %number of time steps could be selected
    %% Proces inputfiles, Copy and make Boundary conditions files depending on the time step.
    try
        % Read mdf-file
        if i == 1
            MDF             = mdf('read',strcat(directory, 'source\',ID1,'.mdf')); %Delft3D mdf file
        else
            MDF             = mdf('read',strcat(directory, 'work\',ID1,'.mdf'));
            NFS             = vs_use(strcat(directory, 'results\trim-', ID1, '_', num2str(i-1),'.dat'),'quiet'); %Use Trim-file previous step
        end
        %Determine the start and stop time from the mdf-file (which is the end time
        %of the previous timestep)
        Time_start      = str2double(MDF.mdf.Data{1,2}(strcmp(MDF.mdf.Data{1,2},'Tstart '),2));
        Time_stop       = str2double(MDF.mdf.Data{1,2}(strcmp(MDF.mdf.Data{1,2},'Tstop  '),2));
        
        % Adjust runtime in mdf-file
        if i == 1
            MDF.mdf.Data{1,2}(strcmp(MDF.mdf.Data{1,2},'Tstart '),2)   = {sprintf('%2.8g',(Time_start))};
            MDF.mdf.Data{1,2}(strcmp(MDF.mdf.Data{1,2},'Tstop  '),2)   = {sprintf('%2.8g',(Time_start + run_timestep))};
        else
            MDF.mdf.Data{1,2}(strcmp(MDF.mdf.Data{1,2},'Tstart '),2)   = {sprintf('%2.8g',(Time_stop))};
            MDF.mdf.Data{1,2}(strcmp(MDF.mdf.Data{1,2},'Tstop  '),2)   = {sprintf('%2.8g',(Time_stop + run_timestep))};
%             if i ==2
%                 a1                          = find(strcmp(MDF.mdf.Data{1,2},'Filic  ') == 1);
%                 MDF.mdf.Data{1,2}(a1,1)     = {strcat('Restid ')}; % Trim-file that gives initial conditions
%             else
                a1                          = find(strcmp(MDF.mdf.Data{1,2},'Restid ') == 1);
%             end
            MDF.mdf.Data{1,2}(a1,2)         = {strcat('trim-', ID1, '_', num2str(i-1))}; % Trim-file that gives initial conditions
            a3                              = find(strcmp(MDF.mor.Data{2,2},'MorStt           ') ==1); % Find morphological start time
            MDF.mor.Data{2,2}(a3,2)         = {'0.0000000e+000      [min]    Spin-up interval from TStart till start of morphological changes'}; % Set morphological start time to 0
            
            copyfile(strcat(directory, 'results\trim-', ID1, '_', num2str(i-1),'.dat'), strcat(directory, 'work\')); % Copy old trimfile to work-directory for initial conditions
            copyfile(strcat(directory, 'results\trim-', ID1, '_', num2str(i-1),'.def'), strcat(directory, 'work\')); % Copy old trimfile to work-directory for initial conditions
        end
        
        % Copy the new data to the work directory (This depends on the
        % model setup
        mdf('write',MDF, strcat(ID1), strcat(directory,'work\'));
        % transport file .tra
        % wind file .wnd
        % config file .xml
        % batch file to run Delft3D .bat
        % restart file  tri-rst
        % observation file .obs
        % roughness file .rgh
        % sediment thickness layer file .sdb
        copyfile(strcat(directory,'source\','***'), strcat(directory, 'work\')); %code to copy file *** have to be filled in for each file       
        
        %% Adjust the boundary files bct and bcc
        INFO = bct_io('READ','*.*.bct'); %original bct file
        for i =1:length(INFO.Table)
            raw_data = INFO.Table(i).Data;
            INFO.Table(i).Data = raw_data(1:14401,:);
            if i==1
                INFO.Table(i).Data=DATA(find(DATA(:,1)>=Time_start & DATA(:,1)<=Time_start+run_timestep),:); % use boundaries for the 
            else
                INFO.Table(i).Data=DATA(find(DATA(:,1)>=Time_stop & DATA(:,1)<=Time_stop+run_timestep),:);
            end
        end
        INFO = bct_io('WRITE','*.*.bct', INFO); 
        
    catch
        disp('error inputfiles');
        return
    end
    
    %% ADD SHOAL MARGING COLLAPSES TO THE DEM (EROSION SCAR & ASSOCIATED DEPOSIT)
    if(i>1)
        try
            NFS             = vs_use(strcat(directory, 'work\','trim-', ID1, '_', num2str(i-1),'.dat'),'quiet'); % Used trim file
            ZZ2=ShoalmarginCollapse3(directory,ID1,i,Ec,Ec2,FSavg,Lsm,GridSize,GridSize2); %Run seperate Matlab script
            %% Change depth in trim -file.
            depth           = vs_get(NFS,'map-sed-series','DPS','quiet');
            ZZ2(isnan(ZZ2))=0;
            depth{end} = -ZZ2;
            NFS=vs_put(NFS,'map-sed-series',{length(depth)},'DPS',{1:length(ZZ2(:,1)) 1:length(ZZ2(1,:))},depth{end}); % Add new Depth file at the end of the TRIM, which will be used as start for the next simulation period
        catch
            disp('error flow-slide')
            return
        end
    end
    
    %% Run simulation
    try
        warning('off')
        run_line = strcat(directory,'work\runmodel.bat');
        cd(strcat(directory, 'work'));
        system(run_line);
        warning('on')
    catch
        disp('error run simulation')
        return
    end
    
    %% Process results into a new folder
    try
        % Copy trim-files and mdf-file to the results-directory
        copyfile(strcat(directory, 'work\trim-', ID1,'.dat'), strcat(directory, 'results\trim-', ID1, '_', num2str(i),'.dat'));
        copyfile(strcat(directory, 'work\trim-', ID1,'.def'), strcat(directory, 'results\trim-', ID1, '_', num2str(i),'.def'));
        copyfile(strcat(directory, 'work\trih-', ID1,'.dat'), strcat(directory, 'results\trih-', ID1, '_', num2str(i),'.dat'));
        copyfile(strcat(directory, 'work\trih-', ID1,'.def'), strcat(directory, 'results\trih-', ID1, '_', num2str(i),'.def'));
        copyfile(strcat(directory, 'work\',      ID1,'.mdf'), strcat(directory, 'results\',      ID1, '_', num2str(i),'.mdf'));
        
    catch
        disp('error process results')
        return
    end
    
end

    