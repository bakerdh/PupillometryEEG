
for subj = 312
tic

    EEGpath = (strcat('/Users/db900/Desktop/Pupildata/EEG/P',num2str(subj),'/'));    
    d_EEG2CSV(EEGpath)
    
toc
end