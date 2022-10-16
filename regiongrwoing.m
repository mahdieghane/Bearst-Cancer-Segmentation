function J=regiongrwoing(I,m,n)
Isizes = size(I); %size of the image
threshI = multithresh(I, 3); %thresholding for three regions
c = impixel(I, m, n); 
currPix = c(1); 
surr = [-1 0; 1 0; 0 -1; 0 1]; %create a mask which represents the four surrounding pixels "4 connection".
mem = zeros(Isizes(1)*Isizes(2), 3); 
mem(1, :) = [m, n, currPix]; 
regSize = 1; %initial size
J = zeros(Isizes(1), Isizes(2)); 
init = 1;
posInList = 1;
k=1; 
%The region growing algorithm.
while(k==1)

    for l=init:posInList 
        for j=1:4 
            m1 = m + surr(j,1);
            n1 = n + surr(j,2);

            check=(m1>=1)&&(n1>=1)&&(m1<=Isizes(1))&&(n1<=Isizes(2)); 

            current = impixel(I, m1, n1);
            currPix = current(1);
            if(check && currPix<=threshI(2) && (J(n1, m1)==0)) %check if it belongs to the thresholding boundary and if not set yet on the image we want to recreate
                posInList = posInList+1;
                mem(posInList, :) = [m1, n1, currPix]; %add the new pixel
                J( n1,m1) = 1;
            end
        end
    end
    if(posInList == init) %when there is no more pixels to add
        k = 0; %make k=0 to close the loop
    else
        init = init+1;
        m = mem(init, 1, :);
        n = mem(init, 2, :);
        k = 1; %keep running the loop
    end
end
end
