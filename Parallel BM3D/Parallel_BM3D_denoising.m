clear all
clc

%This code requires the parallel processing toolbox to increase execution speed
%When not available the parfor loops can be easily
%converted to regular for loops without losing functionality

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMAGE PROCESSING SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampling=13; %bit depth of image, set to 8 for 8 bit images, 13 for 13bit images
             %read the accompanying xml files to be sure of the correct value
             %do not rely on the image bit depth as reported by ImageJ
batch=0;  %Set to 0 to process just one file, set to 1 for processsing multiple files in a folder
stacks=1; %set to 0 for processing single images, set to 1 for stacks
peak=2^sampling; %peak value at a certain sampling bit depth
use_bm3d=1; %Set to 1 to use BM3D denoising, set 0 to use deep convolutional network denoising DnCNN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLAHE SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CLAHE - Contrast Limited Adaptive Histogram Equalization
%Improves image contrast - Slower execution speed
%Do not use for quantitative measurements of fluorescence intensity 
%Useful for linescan images or zstacks
UseCLAHE=0; %set to 1 to enable and 0 to disable
NumTiles=[16,16]; %Number of rectangular contextual regions (tiles) into which adapthisteq divides the image
ClipLimit= 0.01; %Contrast enhancement limit
NBins=peak; %Number of histogram bins used to build a contrast enhancing transformation
Range=('full'); % 'full' to use  the full range of the output class (e.g. [0 255] for uint8) or 'original' to limit the range to [min(I(:)) max(I(:))].




if batch==1
    inpath= uigetdir(pwd, 'Select an input folder');
    images=dir(fullfile(inpath,'*.tif'));  
else
    [singlename, singlepath] = uigetfile('*.tif', 'Select a file folder');
    images=dir(fullfile(singlepath, singlename));
end
    
    outpath=uigetdir(pwd, 'Select an output folder');
    NumberOfImages=length(images); 
    filenames=[images.name];

tic;


if stacks==0 %multiple single images denoised in parallel
    parfor f = 1:NumberOfImages 

        CurrentFileName=getfield(images(f), 'name');
        CurrentFilePath=getfield(images(f), 'folder');
        fileName = strcat(CurrentFilePath,'/',CurrentFileName);
        denoisedfileName= strcat(outpath,'/',CurrentFileName,'_denoised.tif');
        
        disp(sprintf('Processing %s',  fileName));        
        noisy_img = imread(fileName, 1) ; % read in first image
        noisy_img= (double(noisy_img)/peak);
            
            filtered=function_denoise_img(noisy_img,use_bm3d);
            denoised_img=uint16(filtered*peak);
            disp(sprintf('Processed file %s %d/%d', fileName,f,NumberOfImages));
            if UseCLAHE==1
            denoised_img=adapthisteq(denoised_img,'NumTiles', NumTiles, 'ClipLimit', ClipLimit, 'NBins', NBins, 'Range', Range)
            end
            imwrite(denoised_img, denoisedfileName , 'WriteMode' , 'append') ;
               
    end
else %stacks are serially loaded and their frames are denoised in parallel
    for f = 1:NumberOfImages

        CurrentFileName=getfield(images(f), 'name');
        CurrentFilePath=getfield(images(f), 'folder');
        fileName = strcat(CurrentFilePath,'/',CurrentFileName);
        denoisedfileName= strcat(outpath,'/',CurrentFileName,'_denoised.tif');
        
        disp(sprintf('Processing %s',  fileName));
        InfoImage=imfinfo(fileName);
        mImage=InfoImage(1).Width;
        nImage=InfoImage(1).Height;
        frames=length(InfoImage);
        noisy_img=zeros(nImage,mImage,frames,'double'); %creates empty matrix with same dimensions as image to be read
        
        TifLink = Tiff(fileName, 'r'); 
        %frames = numel(imfinfo(fileName)); % return tiff structure, one element per image
        
        %noisy_img = imread(fileName, 1) ; % read in first image
        %noisy_img= (double(noisy_img)/peak);
         
            for ii = 1 : frames         %concatenate each successive tiff to tiff_stack
            TifLink.setDirectory(ii);
            noisy_img(:,:,ii)=TifLink.read();
            %temp_tiff = imread(fileName, ii);
            %temp_tiff= (double(temp_tiff)/peak);
            %noisy_img = cat(3 , noisy_img, temp_tiff);
            end
           
        TifLink.close();
        noisy_img=(noisy_img)/peak;
        
        
        parfor g=1 : frames
            filtered(:,:,g)=function_denoise_img(noisy_img(:,:,g),use_bm3d);
            denoised_img(:,:,g)=uint16(filtered(:,:,g)*peak);
            if UseCLAHE==1
            denoised_img(:, :, g)=adapthisteq(denoised_img(:, :, g),'NumTiles', NumTiles, 'ClipLimit', ClipLimit, 'NBins', NBins, 'Range', Range)
            end
            disp(sprintf('Denoised frame %d/%d of %s' ,g,frames, fileName));
        end
        
        
        for k=1:length(denoised_img(1, 1, :))
            imwrite(denoised_img(:, :, k), denoisedfileName , 'WriteMode' , 'append') ;
        end
        clc
        disp(sprintf('Processed file %s %d/%d', fileName,f,NumberOfImages));

end

    
end
    elapsedtime=round(toc);
    clc
    disp(sprintf('Finished processing all image(s) in %d seconds',elapsedtime));

