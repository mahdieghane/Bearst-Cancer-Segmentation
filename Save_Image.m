function Save_Image(W,L,input_matrix,address)
fid2 = fopen(address,'w');
fprintf(fid2,'P5\n');
fprintf(fid2,'%d %d\n', W+2,L);
fprintf(fid2,'255\n');
% Write the pixel values from the matrix into the file
    y = 1;
    for y =1:L

       for x=1:(W+1)
            if (x==(W+1))
                fwrite(fid2, '\n');
            else
                fwrite(fid2,input_matrix(y,x));
            end
       end
       x = 1;
    end
% Close the file
fclose(fid2);
end