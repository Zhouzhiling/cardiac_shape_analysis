filePath = './dcm_file/';
outputPath = './dcm2png/';
imagetype = 'png';

listing = dir(filePath);
numofCase = length(listing);

for k = 3:numofCase
    caseName = listing(k).name;
    outPath = strcat(outputPath,caseName);
    mkdir(outPath);
    
    % get the sequence number of corresponding txt file
    txtPath = strcat('./txt/',caseName,'/contours-manual/IRCCI-expert/');
    listingTXT = dir(txtPath);
    numofTXT = length(listingTXT);
    seqArr = zeros(1,numofTXT-2);
    
    % record the sequence array
    for i = 3 : numofTXT
        TXTName = listingTXT(i).name;
        seq = str2num(TXTName(9:12));
        seqArr(i-2) = seq;
    end
    seqArr = unique(seqArr);
    
    dcmPath = strcat(filePath,caseName);
    listingD = dir(dcmPath);
    numofDCM = length(listingD);
    for i = 3 : numofDCM
        DCMName = listingD(i).name;
        seq = str2num(DCMName(9:12));
        if(find(seqArr == seq)) % only output the png file that has the corresponding txt contour file
            dcmFile = strcat(dcmPath,'/',DCMName);
            info = dicominfo(dcmFile);
            dcmfile = dicomread(info);
            figure('visible','off')
            imshow(dcmfile,[]);
            dcmImagei = uint8(255 * mat2gray(dcmfile)); %Convert to uint8 format
            filename = strcat(outPath, '/', DCMName(1:end-4),'.png');
            imwrite(dcmImagei,filename, imagetype);% Save Image to specified image type
            display(filename);
        end
    end
end