function project_main()
    close all;
    clear all;
    % Loaded image resized into a square pixels
    n = 400;

    folder = 'tools/';

    direc = dir(strcat(folder,'*.bmp'));

    % Set time parameters
    % tStep = 0.1;
    iter = 2000;

    % Set scaling parameters for energy terms
    lambda= 0.25;

    % nContours: number of contours in each row column
    nContours = 1;
    radius = n/nContours/4;

    % Noise trigger: if 'true' noise added
    % Noise types 1:gaussion, 2:salt & pepper, 3:localvar, 4:poisson, 5:speckle
    noiseTrigger = true;
    noiseAll = {'gaussian'; 'salt & pepper'; 'poisson'; 'speckle'};
    noiseType = noiseAll{2};

    % Create two matrices where each point has the value of the coordinate
    % split into corresponding positions in X and Y matrices
    X = zeros(n,n);
    Y = zeros(n,n);
    for i = 1 : n
        for j = 1 : n
            X(i,j) = j;
            Y(i,j) = i;
        end
    end

    % Generate contours with equally spaced centers
    k=1;
    psiInit(1,:,:)=zeros(n,n);
    for i = 1:nContours
        for j = 1:nContours
            f1=1+(i-1)*2;
            f2=1+(j-1)*2;
            centerX = n/nContours/2*f1;
            centerY = n/nContours/2*f2;
            psiInit(k,:,:) = sqrt((X-centerX).^2 + (Y-centerY).^2)-radius;
            k=k+1;
        end
    end


    for cur_noise=1:length(noiseAll)+1
        if cur_noise > length(noiseAll)
            noiseTrigger=false;
        else
            noiseType=noiseAll{cur_noise};
        end


        % Load image: Change filename with director to image file
        filename = direc(1).name;
        image = imread(strcat(folder,filename));

        % Convert to greyscale if RGB
        if size(image,3)==3
            image = rgb2gray(image);
        end


        % Resize image and nomalize
        image = imresize(image, [n, n]);
        image = double(image);
        image = image./max(image(:));

        % Average image intensity of max and min
        % image_avg = sum(image(:))/(n*n);
        image_avg = (max(image(:))+min(image(:)))/2;

        % Pixels having intensity higher than average
        radius_dyn = sum(image(:)>image_avg);

        % radius: radius of each contour
        % radius=radius_dyn;

        cleanI = image;


        % Add noise if specified and normalize as chosen in specifications
        if noiseTrigger
            image = imnoise(image, noiseType);
        end


        % Find the minimum of each point to get overall initial contour
        psi = squeeze(min(psiInit,[],1));

        % Compute segmentation with Chan_Vese
        I = image;


        for i = 1:iter
        %     Average intensity inside contour
            u = sum(sum(I.*(psi < 0)))/sum(psi(:) < 0);

        %     Average intensity outside contour
            v = sum(sum(I.*(psi > 0)))/sum(psi(:) > 0);   

        %     Obtain tStep from CFL condition
            tStep = min( 0.5/(max(I(:).^2)*(1-lambda)), 1/(lambda*sqrt(2)) );

        %     Determine the psi at the next time step using the update equation
            psi = psi +  tStep.*((1-lambda)*((I-u).^2.*normal(psi,1)...
                - (I-v).^2.*normal(psi,-1)) + lambda.*kappaNormal(psi)) ;

            acc(i) = sum(sum((psi>0) == (cleanI>0.5)))/(n*n);
        end
        
        figure(cur_noise);
        plot(1:iter,acc);
        avg_acc(cur_noise) = sum(acc)/length(acc);
    end
    
    disp(avg_acc)
end

% Compute the Level Set formulation of normal using upwind entropy scheme 
% and edge preservation
function norm = normal(psi, sign)
% Preserve the edges so identical values at extension of image edge
    DxF = psi([2:end end],:)-psi;
    DxB = psi - psi([1 1:end-1],:);
    DyF = psi(:,[2:end end])-psi;
    DyB = psi - psi(:,[1 ,1:end-1]);

% Level set formulation of normal
    if sign == 1
        norm = ((max(DxF,0)).^2 + (min(DxB,0)).^2 + (max(DyF,0)).^2 ...
            + (min(DyB,0)).^2).^(1/2);
    elseif sign == -1
        norm = ((min(DxF,0)).^2 + (max(DxB,0)).^2 + (min(DyF,0)).^2 ...
            + (max(DyB,0)).^2).^(1/2);
    end
end

% Compute the level set formulation of kappa * Normal using a central
% difference scheme and edge preservation
function kn = kappaNormal(psi)
% Preserve the edges so identical values at extension of image edge
    psix = ( psi([2:end end],:)-psi([1 1:end-1],:) )./2;
    psiy = ( psi(:,[2:end end])-psi(:,[1 ,1:end-1]) )./2;
    psixy = ( psix(:,[2:end end])-psix(:,[1 ,1:end-1]) )./2;
    psixx = ( psi([2:end end],:)- 2.*psi + psi([1 1:end-1],:) );
    psiyy = ( psi(:,[2:end end])- 2.*psi + psi(:,[1 ,1:end-1]) );
    
% Level set formulation of kappa*normal
    kn = (psixx.*(psiy.^2) - 2.*psix.*psiy.*psixy + psiyy.*(psix.^2))...
        ./(max(eps, ((psix.^2)+(psiy.^2))));
end