function d_EEG2CSV(EEGpath)

% function to convert ANT EEProbe files to (zipped) csv format
% uses an EEGlab plugin to read the files, so EEGlab must be installed and on the path
% input is the directory containing raw EEG files - it is critical that the entire path is given, the EEGlab functions do not understand the shortcut '~/' to mean the home directory
% each .cnt file will be saved as a gzipped csv file in the same directory
% these files are up to 5x larger than the originals!
% the columns in the csv file are Time (ms), Trigger codes, followed by all electrodes included in the montage
% rows are sample points
% ideally these files should be accompanied by meta-data describing the experiment
% DHB 7/9/18

d = dir(EEGpath);       % index the directory containing the EEG files
counter = 0;
for n = 1:length(d)
    temp = d(n).name;
    if length(temp)>3
        if temp(end-2:end)=='cnt'       % find each file in cnt format
            counter = counter + 1;
            namelist{counter} = temp;   % create a list of cnt filenames
        end
    end
end

for c = 1:counter                       % loop through all cnt files
    fname = namelist{c};
    EEG = pop_loadeep_v4(strcat(EEGpath,fname));    % load the EEG data
    %     save(strcat(EEGpath,fname(1:end-3),'mat'),'EEG');     % optionally save the EEG structure in .mat format

    eventcounter = 0;
    alllatency = [];
    alltypes = [];    
    for n = 1:length(EEG.event)                     % convert trigger information to a useable format
        if length(EEG.event(n).type)<4
            eventcounter = eventcounter + 1;
            alllatency(eventcounter) = EEG.event(n).latency;
            alltypes(eventcounter) = str2num(EEG.event(n).type);
        end
    end
    outputfile = strcat(EEGpath,fname(1:end-3),'csv');  % generate a filename for the output based on the original filename
    fid = fopen(outputfile,'w');                        % open the output file
    electrodenames = '';                                % initialise a string to contain the names of all the electrodes
    formatvector = '%2.0f, %2.0f';                      % initialise a string determining the format of the numbers we will output
    for ch = 1:EEG.nbchan           % inside this loop, add values to the above strings, depending on how many channels there are
        electrodenames = strcat(electrodenames,EEG.chanlocs(ch).labels,', ');
        formatvector = strcat(formatvector,', %2.6f');
    end
    fprintf(fid,strcat('Time, Trigger, ',electrodenames(1:end-1),'\n'));    % output the column names to the text file
    for t = 1:EEG.pnts                              % loop through all samples
        trcode = 0;                                 % default state is no trigger
        if ismember(EEG.times(t),alllatency)        % check if a trigger occurred at this time point
            i = find(EEG.times(t)==alllatency);
            trcode = alltypes(i);                   % if it did, output the trigger code at this time point
        end
        fprintf(fid,strcat(formatvector,'\n'),[EEG.times(t),trcode,EEG.data(:,t)']);    % export this row of data to the text file
    end
    fclose(fid);                    % close the file
    gzip(outputfile);               % gzip the file to save space (usually about 60% smaller)
    delete(outputfile);             % delete the original csv file
    clear EEG alllatency alltypes;  % clean up variables
end

end