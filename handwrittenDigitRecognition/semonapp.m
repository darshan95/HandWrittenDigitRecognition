% clc;
clear;
close all;
% Num : 0-11
%      5
%      0
%      4
%      1
%      9
%      2
%      1
%      3
%      1
%      4
%      3
% Num : 18 -> 6 (bad)
% Num : 15 -> 7
% Num : 17,31 -> 8
% Num : 33 -> 4 (closed)

tot_num_imgs =700 ;
offset = 1;
ans=0;

corr=zeros(11);
graph=zeros(11);

[img_tot lab_tot] = readMNIST('train-images.idx3-ubyte','train-labels.idx1-ubyte', tot_num_imgs,offset);
%lab_t,ot

actual=zeros(1,11);
found=zeros(1,11);
r=0;

for img_num=1:tot_num_imgs
    img = img_tot(:,:,img_num);
    lab = lab_tot(img_num);
    [m n] = size(img);
    actual(1,lab+1)=actual(1,lab+1)+1;
    %     figure,
    %     imshow(img);
    %BW3 = bwmorph(img,'remove');
    %figure,
    %imshow(BW2);
    %      BW3 = bwmorph(img,'skel',Inf);
    %      se = strel('disk',2);
    
    %      BW3=imdilate(BW3,se);
    %      BW3=imerode(BW3,se);
    %figure,
    %imshow(BW3);
    %     img=BW3;
      
    % Seems to be making things worse. Can chuck this out.
%     if lab==8
%         figure,
%         subpolot(2,4,1)
%         imshow(img);xlabel('original');
%     end
    
    %     %Select disk-SE of size > width of line
%          se = strel('disk',1);
    %     %Perform Closing of image.
%          img = imclose(img,se);
    
    % Binary and Complement.
  
    img = im2bw(img);

    img_comp = imcomplement(img);
 
    
    %Filling the holes.
    img_fill = img_comp;
    con_comps = bwconncomp(img_fill);
    numPixels = cellfun(@numel,con_comps.PixelIdxList);
    [biggest,idx] = max(numPixels);
    for i=1:con_comps.NumObjects
        if i ~= idx
            img_fill(con_comps.PixelIdxList{i}) = 0;
        end
    end
    
    %Complement again
    img_fill_comp = imcomplement(img_fill);
    

    % Taking difference of images to locate blobs.
    blob_img = img;
    for i=1:m
        for j=1:n
            blob_img(i,j) = img_fill_comp(i,j) - img(i,j);
            if blob_img(i,j) < 0
                blob_img(i,j) = 0;
            end
        end
    end
    
    % Additional : Makes the blob cover actual blob size so finding stems
    % is easier.
    
%     se = strel('square',1);
%     blob_img = imdilate(blob_img,se);
    
    % Find the number of blobs.
    blob_cc = bwconncomp(blob_img);
    %removing additional noise
    thresh = 2;
    for i=1:blob_cc.NumObjects
        if (bwarea(blob_img(blob_cc.PixelIdxList{i})) < thresh)
            blob_img(blob_cc.PixelIdxList{i}) = 0;
        end
    end
    
    stem_img = img;
    
    % Group 1 (digit-8)
    if blob_cc.NumObjects == 2
        digit = 8;
        
        % Group 2 (digit-4,6,9,0)
    elseif blob_cc.NumObjects == 1
        %     elseif blob_cc.NumObjects == 1
        for i=1:m
            for j=1:n
                stem_img(i,j) = img_fill_comp(i,j) - blob_img(i,j);
                if stem_img(i,j) < 0
                    stem_img(i,j) = 0;
                else
                    
                end
            end
        end
        % Additional : Removing the noise.
     stem_cc = bwconncomp(stem_img);
        
        %again because image has noise removed now.
        thres=4;
        if stem_cc.NumObjects == 1 % Num is 6 or 9 (locate stem-blob relative position)
            blob_cent = regionprops(blob_cc,'centroid');
            stem_cent = regionprops(stem_cc,'centroid');
            if  (blob_cent.Centroid(2)-stem_cent.Centroid(2))^2 + (blob_cent.Centroid(1)-stem_cent.Centroid(1))^2 < thres
                digit=0;
            elseif blob_cent.Centroid(2) < stem_cent.Centroid(2) % because (0,0) is top-left.
                digit = 9;
            else
                digit = 6;
            end
        else % Num is 4
            digit = 4;
        end
        
        % Group 3 (digit-4,1,2,3,5,7)
    else
        digit = 412357;
        %          removing additional noise
        img_cc = bwconncomp(img);
        thresh = 2;
        for i=1:img_cc.NumObjects
            if (bwarea(img(img_cc.PixelIdxList{i})) < thresh)
                img(img_cc.PixelIdxList{i}) = 0;
            end
        end
        br_fl = 0;
        for j=1:n
            for i=1:m/2
                if img(i,j) == 1 && img(i,j+1)==1
                    top_left = [i,j];
                    br_fl = 1;
                    break;
                end
            end
            if br_fl == 1
                break;
            end
        end
        br_fl = 0;
        for j=1:n
            for k=1:m/4
                k=m-k;
                if img(k,j) == 1 && img(k,j+1)==1
                    bot_left = [k ,j];
                    br_fl = 1;
                    break
                end
            end
            if br_fl == 1
                break
            end
        end
        %Drawing a line from top left position to bottom left
        K=img;
        rpts = linspace(top_left(1),bot_left(1),1000);   %# A set of row points for the line
        cpts = linspace(top_left(2),bot_left(2),1000);   %# A set of column points for the line
        index = sub2ind([m n],round(rpts),round(cpts));  %# Compute a linear index
        K(index) = 1;
%         figure,imshow(K);
        
        %  K = im2bw(K);
        
        
        %          subplot(5,5,img_num);
        %         imshow(K);
        %finding the Blob
        % Binary and Complement.
        
        
        % To connect disconnected image.
        %se = strel('disk',1);
        %  img = imclose(img,se);
        %se1 = strel('square',1);
        %img = imopen(img,se);
        K_comp = imcomplement(K);
        %                 se = strel('disk',1);
        %      K_comp = imclose(K_comp,se);
        
        %Filling the holes.
        K_fill = K_comp;
        con_comps_k = bwconncomp(K_fill);
        numPixels = cellfun(@numel,con_comps_k.PixelIdxList);
        [biggest,idx] = max(numPixels);
        for i=1:con_comps_k.NumObjects
            if i ~= idx
                K_fill(con_comps_k.PixelIdxList{i}) = 0;
            end
        end
        %Complement again
        K_fill_comp = imcomplement(K_fill);
        
        
        % Seems to be making things worse. Can chuck this out.
        
        %     %Select disk-SE of size > width of line
        %     se = strel('disk',2);
        %     %Perform Closing of image.
        %     img_closed = imclose(img,se);
        
        % Taking difference of images to locate blobs.
        blob_K = K;
        for i=1:m
            for j=1:n
                blob_K(i,j) = K_fill_comp(i,j) - K(i,j);
                if blob_K(i,j) < 0
                    blob_K(i,j) = 0;
                end
            end
        end
        
        % Additional : Makes the blob cover actual blob size so finding stems
        % is easier.
        
        %   se = strel('disk',3);
        %  blob_K = imdilate(blob_K,se);
        
        % Find the number of blobs.
        blob_kk = bwconncomp(blob_K);
        %removing additional noise
        thresh = 10;
        for i=1:blob_kk.NumObjects
            if (bwarea(blob_K(blob_kk.PixelIdxList{i})) < thresh)
                blob_K(blob_kk.PixelIdxList{i}) = 0;
            end
        end
        
        if blob_kk.NumObjects == 0
            L=img;
            br_fl = 0;
            for j=0:n-1
                j=n-j;
                for k=0:m/4
                    k=m-k;
                    if img(k,j) == 1
                        bot_right = [k ,j];
                        br_fl = 1;
                        break
                    end
                end
                if br_fl == 1
                    break
                end
            end
            rpts = linspace(top_left(1),bot_right(1),1000);   %# A set of row points for the line
            cpts = linspace(top_left(2),bot_right(2),1000);   %# A set of column points for the line
            index = sub2ind([m n],round(rpts),round(cpts));  %# Compute a linear index
            L(index) = 1;
            
            
            % To connect disconnected image.
            %se = strel('disk',1);
            %  img = imclose(img,se);
            %se1 = strel('square',1);
            %img = imopen(img,se);
            L_comp = imcomplement(L);
            
            %Filling the holes.
            L_fill = L_comp;
            con_comps = bwconncomp(L_fill);
            numPixels = cellfun(@numel,con_comps.PixelIdxList);
            [biggest,idx] = max(numPixels);
            for i=1:con_comps.NumObjects
                if i ~= idx
                    L_fill(con_comps.PixelIdxList{i}) = 0;
                end
            end
            
            %Complement again
            L_fill_comp = imcomplement(L_fill);
            
            % Seems to be making things worse. Can chuck this out.
            
            %     %Select disk-SE of size > width of line
            %     se = strel('disk',2);
            %     %Perform Closing of image.
            %     img_closed = imclose(img,se);
            
            % Taking difference of images to locate blobs.
            blob_L = L;
            for i=1:m
                for j=1:n
                    blob_L(i,j) = L_fill_comp(i,j) - L(i,j);
                    if blob_L(i,j) < 0
                        blob_L(i,j) = 0;
                    end
                end
            end
            
            % Additional : Makes the blob cover actual blob size so finding stems
            % is easier.
            
            %se = strel('disk',3);
            %blob_L = imdilate(blob_L,se);
            
            % Find the number of blobs.
            blob_ll = bwconncomp(blob_L);
            %removing additional noise
            thresh = 2;
            for i=1:blob_ll.NumObjects
                if (bwarea(blob_L(blob_ll.PixelIdxList{i})) < thresh)
                    blob_L(blob_ll.PixelIdxList{i}) = 1;
                end
            end
            if blob_ll.NumObjects == 0
                digit=1;
            else
                digit=4;
            end
            
        else
            %se = strel('disk',3);
            %blob_K = imdilate(blob_K,se);
            stem_K=K;
            for i=1:m
                for j=1:n
                    stem_K(i,j) = K_fill_comp(i,j) - blob_K(i,j);
                    if stem_K(i,j) < 0
                        stem_K(i,j) = 0;
                    end
                end
            end
            stem_kk = bwconncomp(stem_K); %again because image has noise removed now.
            if stem_kk.NumObjects == 1
                blob_cent_k = regionprops(blob_kk,'centroid');
                stem_cent_k = regionprops(stem_kk,'centroid');
                [s11,s12]=size(blob_cent_k);
                [s21,s22]=size(stem_cent_k);
%                 if  ( blob_cent_k.Centroid(1,2) < stem_cent_k.Centroid(1,2) ) % because (0,0) is top-left.
%                     digit = 2;
%                 else
                    digit=5;
%                 end
            else
                br_fl = 0;
                for j=0:n-1
                    j=n-j;
                    for i=1:m/2
                        if stem_K(i,j) == 1  && img(i,j-1)==1
                            top_right = [i,j];
                            br_fl = 1;
                            break;
                        end
                    end
                    if br_fl == 1
                        break;
                    end
                end
                br_fl = 0;
                for j=0:n-1
                    j=n-j;
                    for k=0:m/4
                        k=m-k;
                        if stem_K(k,j) == 1
                            bot_right = [k ,j];
                            br_fl = 1;
                            break
                        end
                    end
                    if br_fl == 1
                        break
                    end
                end
                rpts = linspace(top_right(1),bot_right(1),1000);   %# A set of row points for the line
                cpts = linspace(top_right(2),bot_right(2),1000);   %# A set of column points for the line
                index = sub2ind([m n],round(rpts),round(cpts));  %# Compute a linear index
                stem_K(index) = 1;
                
                stem_K_comp = imcomplement(stem_K);
                
                %Filling the holes.
                stem_K_fill = stem_K_comp;
                con_comps = bwconncomp(stem_K_fill);
                numPixels = cellfun(@numel,con_comps.PixelIdxList);
                [biggest,idx] = max(numPixels);
                for i=1:con_comps.NumObjects
                    if i ~= idx
                        stem_K_fill(con_comps.PixelIdxList{i}) = 0;
                    end
                end
                
                %Complement again
                stem_K_fill_comp = imcomplement(stem_K_fill);
                
                % Seems to be making things worse. Can chuck this out.
                
                %     %Select disk-SE of size > width of lin              %     se = strel('disk',2);
                %     %Perform Closing of image.
                %     img_closed = imclose(img,se);
                
                % Taking difference of images to locate blobs.
                stem_K_blob = stem_K;
                for i=1:m
                    for j=1:n
                        stem_K_blob(i,j) = stem_K_fill_comp(i,j) - stem_K(i,j);
                        if stem_K_blob(i,j) < 0
                            stem_K_blob(i,j) = 0;
                        end
                    end
                end
                
                % Additional : Makes the blob cover actual blob size so finding stems
                % is easier.
                
                se = strel('disk',3);
                stem_K_blob = imdilate(stem_K_blob,se);
                
                % Find the number of blobs.
                blob_kk = bwconncomp(stem_K_blob);
                %removing additional noise
                thresh = 2;
                for i=1:blob_kk.NumObjects
                    if (bwarea(stem_K_blob(blob_kk.PixelIdxList{i})) < thresh)
                        stem_K_blob(blob_kk.PixelIdxList{i}) = 1;
                    end
                end
                if blob_kk.NumObjects == 0
                    digit=7;
                else
                    digit=3;
                end
            end
        end
    end
    if lab==2
     figure;
    
        subplot(2,4,2);
        imshow(img);xlabel('Dilated');
        subplot(2,4,3);
        imshow(img_comp);xlabel('complement of original');
        subplot(2,4,4);
        imshow(img_fill);xlabel('complement filled');
        subplot(2,4,5);
        imshow(img_fill_comp);xlabel('again complement');
        subplot(2,4,6);
        imshow(blob_img);xlabel('blobs');
        subplot(2,4,7);
        imshow(stem_img);xlabel('stems');
    end
    
    %  break
    %
    % digit
    if digit==lab
        ans=ans+1;
        found(1,digit+1)=found(1,digit+1)+1;
    end
    
    graph(lab+1,digit+1)=graph(lab+1,digit+1)+1;
    
    %  subplot(5,5,img_num);
    %    imshow(img);
    %  xlabel(lab);
    %   ylabel(digit);
    %         if(digit==412357)
    %         hold on;
    %         imshow(img); plot([top_left(2),bot_left(2)],[top_left(1),bot_left(1)],'Color','r','LineWidth',2);
    %         hold off;
    %     end
    
    %     break
    
    % Closing with disk size greater than width of lines.
    
    
    %     bw = bwconncomp(img);
    %
    %     subplot(4,4,img_num);
    %     imshow(img);
    %     xlabel(lab);
    %     num_blobs = bw.NumObjects-1;
    %     if (num_blobs==2)
    %         digit = 8;
    %     elseif (num_blobs==1)
    %         i1=imcomplement(img);
    %         i2 = imfill(i1,'holes');
    %         i3 = imcomplement(i2);
    %         se=strel('disk',3);
    %         i4 = imclose(i3,se);
    %         i5 = i4 - i3;
    %         blob_img = imcomplement(i5);
    %         %need to use bwboundaries to find the number of stems and its
    %         %position
    %         digit=9;
    %     elseif (num_blobs==0)
    %         digit = 0;
    %     end
    %     disp(num_blobs);
    %     ylabel(digit);
end
graph
actual
found
per=zeros(1,11);
for i=1:10
    per(i)=(found(1,i)*100.0)/actual(1,i);
end
per
% %r
% %ans
 (ans*100.0)/tot_num_imgs