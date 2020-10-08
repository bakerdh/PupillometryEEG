function d_EEGheader

% function to create a generic meta-data file for an EEG experiment
% contains electrode positions and trigger explanations
% can contain other useful study information
% DHB 30/9/18

load Waveguard.mat;
    outputfile = 'headerfile.csv';  % generate a filename for the output based on the original filename
    fid = fopen(outputfile,'w');                        % open the output file


    detailvect1{1} = 'Binocular flicker pupillometry and EEG';
    detailvect2{1} = '2019';
    detailvect1{2} = 'Sample rate';
    detailvect2{2} = '1000Hz';    
    detailvect1{3} = 'Stimulus duration';
    detailvect2{3} = '12s';
    detailvect1{4} = 'Target frequency';
    detailvect2{4} = '2Hz';
    detailvect1{5} = 'Mask frequency';
    detailvect2{5} = '1.6Hz';
    detailvect1{6} = 'Repetitions per condition';
    detailvect2{6} = '3';
    detailvect1{7} = 'Trials per block';
    detailvect2{7} = '60';
    detailvect1{8} = 'Total participants';
    detailvect2{8} = '30';
    detailvect1{9} = 'Number of electrodes';
    detailvect2{9} = '64';
    detailvect1{10} = 'Montage layout';
    detailvect2{10} = '5 percent';
    detailvect1{11} = 'EEG system';
    detailvect2{11} = 'ANT Neuroscan';
    detailvect1{12} = 'Original file format';
    detailvect2{12} = 'ANT EEprobe .cnt files';
    
    triggervect = 1;
    triggerinfo{1} = 'Binary triggers require additional files for interpretation';

               
    plist = 1:30;
    fprintf(fid,'Details1,Details2,Trigger,Description,ParticipantIDs,Electrode,X_position,Y_position,OutlineX,OutlineY,NoseX,NoseY,LearX,LearY,RearX,RearY\n');
 
 for n = 1:101
     
     if n<=length(detailvect1)
        detailstring1 = detailvect1{n};
     else
         detailstring1 = '';
     end
     if n<=length(detailvect2)
        detailstring2 = detailvect2{n};
     else
         detailstring2 = '';
     end     
     if n<=length(triggervect)
        t = num2str(triggervect(n));
     else
         t = '';
     end   
     if n<=length(triggerinfo)
        triggerstring = triggerinfo{n};
     else
         triggerstring = '';
     end     
     if n<=length(plist)
         s = num2str(300+plist(n));
         while length(s)<2
             s = strcat('0',s);
         end
        pID = strcat('P',s);
     else
         pID = '';
     end    
     
     if n<=length(lay.label)
       electrodelabel = lay.label{n};
       x = num2str(lay.pos(n,1),'%2.4f');
       y = num2str(lay.pos(n,2),'%2.4f');
     else
       electrodelabel = '';
       x = '';
       y = '';
     end
     
         xhead = num2str(lay.outline{1,1}(n,1),'%2.4f');
         yhead = num2str(lay.outline{1,1}(n,2),'%2.4f');
         
     if n<=length(lay.outline{1,2})
         xnose = num2str(lay.outline{1,2}(n,1),'%2.4f');
         ynose = num2str(lay.outline{1,2}(n,2),'%2.4f');
     else
       xnose = '';
       ynose = '';
     end
     if n<=length(lay.outline{1,4})
         xLear = num2str(lay.outline{1,4}(n,1),'%2.4f');
         yLear = num2str(lay.outline{1,4}(n,2),'%2.4f');
     else
       xLear = '';
       yLear = '';
     end
     if n<=length(lay.outline{1,3})
         xRear = num2str(lay.outline{1,3}(n,1),'%2.4f');
         yRear = num2str(lay.outline{1,3}(n,2),'%2.4f');
     else
       xRear = '';
       yRear = '';
     end
     
fprintf(fid,strcat(detailstring1,',',detailstring2,',',t,',',triggerstring,',',pID,',',electrodelabel,',',x,',',y,',',xhead,',',yhead,',',xnose,',',ynose,',',xLear,',',yLear,',',xRear,',',yRear,'\n'));    % export this row of data to the text file

 end
 fclose(fid); 
 
end