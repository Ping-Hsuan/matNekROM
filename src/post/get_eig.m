gu = dlmread("../ops/gu");
gu = reshape(gu,sqrt(length(gu)),sqrt(length(gu)));

[V, D] = eig(gu);

eig = sort(diag(D),'descend');
fileID = fopen('eig.txt','w');
fprintf(fileID,'%24.16f\n',eig);
fclose(fileID);
