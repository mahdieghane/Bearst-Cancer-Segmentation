
close all;
clear all;

Test_Image_1=imread('Test_Image_20.png');
Mask_1=imread('mask_20.png');

[m,n]=size(Test_Image_1);
%m=128;n=128;
%Test_Image_1=zeros(m,n);
%Test_Image_1=Test_Image(1:128,1:128,1);
%Mask_1=zeros(m,n);
%Mask_1=Mask(1:128,1:128,1);

Test_Image_1=imgaussfilt(Test_Image_1,0.3);
med_img = medfilt2(Test_Image_1);

hist_img = histeq(med_img);  
%imshow(hist_img,[])

%% opening and erosion
se = strel('square',9);
Ie = imerode(hist_img,se);
Io = imopen(Ie,se);
Iobr = imreconstruct(Io,hist_img);
%figure();imshow(Iobr)
%title('Opening-by-Reconstruction')
%% closing and dialation
se1 = strel('diamond',1);
Ioc = imclose(Io,se1);
%imshow(Ioc)
Iobrd = imdilate(Iobr,se1);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
%figure();imshow(Iobrcbr)
%title('Closing by Reconstruction')
%% Threshold
%level = graythresh(hist_img);
%BW = imbinarize(hist_img,level);
BB=Iobrcbr<40;
img_binary = imbinarize(Iobrcbr);
%figure(); imshow(img_binary)

%% More dilation and erosion
se1 = strel('disk', 3);  
se2 = strel('square', 3);
dilated_img = imdilate(BB, se1); 
erode_img = imerode(dilated_img, se2);
%figure(); imshow(erode_img)

%% Dilation and Erosion on Complete image
comp_img_ = imcomplement(erode_img); 
%figure(); imshow(comp_img_)
se3 = strel('square', 3);  
se4 = strel('diamond', 3);
dilated_img_1 = imdilate(comp_img_, se3); 
erode_img_1 = imerode(dilated_img_1, se4);
%figure(); imshow(erode_img_1)

comp_img_c = imcomplement(erode_img_1); 
clear_img = imclearborder(comp_img_c, 4); 
%figure(); imshow(clear_img)

%%

bwar_img = bwareafilt(clear_img, 1); 
%figure(); imshow(bwar_img)
%%
out = imfill(bwar_img,'holes');
%figure(); imshow(out)
CH = bwconvhull(out,'objects');
%imshow(CH,[])
%% intersection
edgeD=edge(hist_img,'canny',0.2);
%imshow(edgeD,[])

Div=zeros(m,n);
for i=1:m
    for j=1:n
        if CH(i,j)==1 && edgeD(i,j)==1
            Div(i,j)=255;
        end
    end
end
%figure();imshow(Div,[])
%%
% [row,column]=find(Div~=0);
% res=regiongrwoing(hist_img,row(80,1),column(80,1));
% figure();imshow(res,[]);
%% More dialation

se4=strel('square',3);
se5=strel('disk',1);
se6=strel('disk',3);

Dialated_Div=imdilate(Div,se4);
%figure();imshow(Dialated_Div,[])
Eroded_Div=imopen(Dialated_Div,se5);
%figure();imshow(Eroded_Div,[])
Dialated_Div1=imdilate(Eroded_Div,se4);
Dialated_Div2=imdilate(Dialated_Div1,se4);
Dialated_close = imclose(Dialated_Div2,se6);
%figure();imshow(Dialated_close,[])

Filled=imfill(Dialated_close,4);
%figure();imshow(Filled,[]);
Filled1=imfill(Filled,8);
%figure();imshow(Filled1,[]);
Convexed1 = bwconvhull(Filled1);
figure();imshow(Convexed1);
title('Result of edge detection and dilation erosion')
%%Accuracy
Accuracy = jaccard(logical(Mask_1), Convexed1);
%% Region Growing

region_center=regionprops(Convexed1,'Centroid');
C=struct2cell(region_center);
A=cell2mat(C);
Out=regiongrwoing(hist_img,round(A(1,1)),round(A(1,2)));
Output=Out.*double(Convexed1);
%Ac=jaccard(logical(Mask_1), logical(Output));
F=zeros(m,n);
for i=1:m
    for j=1:n
    if Output(i,j)==1 || Div(i,j)==1
        F(i,j)=255;
    end
    end
end
%figure();imshow(F,[]);
Filled_i=imfill(F,8);
%figure();imshow(Filled_i,[]);
%% final result
se8=strel('disk',1);
se9=strel('disk',3);

F_eroded=imerode(Filled_i,se8);
F_dilated=imdilate(F_eroded,se9);
figure();imshow(F_dilated);
title('result of Region Grwoing')
Ac_growing_region=jaccard(logical(Mask_1), logical(F_dilated));

%% gaussian filtering (noise removal)
img_gaussian = imgaussfilt(Test_Image_1, 0.82);

%% watershed

I1 = imtophat(img_gaussian, strel('disk', 70));
I2 = imadjust(I1);
level = graythresh(I2);
BW = im2bw(I2, level);
C=BW;
D = -bwdist(C);
D(C) = - Inf;
L = watershed(D);
water_hist = histeq(L);
waterhist=imclearborder(water_hist,8);
%%
X=double(waterhist).*double(Convexed1);
%figure();imshow(X,[]);
%%
threshold=0.3;
eshtrak=zeros(m,n);
t=0;
r=0;
s=0;
mini=min(min(waterhist));
for i=1:m
    for j=1:n
        if X(i,j)~=mini
            s=find(X==X(i,j));
            [t,r]=find(waterhist==X(i,j));
            if ((length(s))/length(t))>threshold
               for k=1:length(t)
                   eshtrak(t(k),r(k))=1;
               end
            end
            t=0;
            r=0;
        end
    end
end
eshtrak=255*eshtrak;
%figure();imshow(eshtrak,[])
%%
clear_w=imclearborder(eshtrak,4);
dilated_eshtrak=imdilate(clear_w,se9);
open_d=imopen(dilated_eshtrak,se8);
%figure();imshow(dilated_eshtrak,[]);
%figure();imshow(open_d,[]);
fgm = bwareafilt(logical(open_d),1);
%figure();imshow(fgm,[]);
Ac_watershed=jaccard(logical(Mask_1), fgm);
%%
Final_im=zeros(m,n);
for i=1:m
    for j=1:n
        if fgm(i,j)==1 || F_dilated(i,j)==255
            Final_im(i,j)=255;
        end

    end
end
figure();imshow(Final_im,[])
title('Final Result')
Final_accuracy=jaccard(logical(Mask_1), logical(Final_im));
%%
RF=[Ac_growing_region,Ac_watershed,Accuracy,Final_accuracy];
Max_accuracy=max(RF);
if Max_accuracy == Ac_growing_region
        method='region-growing';
elseif Max_accuracy==  Ac_watershed
        method='watershed';
elseif Max_accuracy==    Accuracy
        method='Edge_detection,Erosion_dilation Accuracy';
else  
        method='Union of region_growing and watershed';
end
str=strcat('Maximum accuracy is for method :',method);
Accuracies=[Accuracy;Ac_growing_region;Ac_watershed;Final_accuracy;Max_accuracy]*100;
Methods={'Edge_detection & Erosion_dilation Accuracy';'region-growing Accuracy' ...
          ;'watershed Accuracy';'Union of region_growing and watershed Accuracy';str};
T = table(categorical(Methods),Accuracies,'VariableNames',{'Methods','Accuracies'})
%% Save the result
%Save_Image(m,n,Final_im,'Final_result.pgm');
