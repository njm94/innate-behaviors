function [CM,P] = OIA_SA_seedpixel(I,x,y)
% Matthieu Vanni - Sept 27, 2013
% function [CM,P] = OIA_SA_seedpixel(I,x,y)
% create the matrix of correlation between the seed pixel and every pixel
% of I
% I: matrix of data
% x,y: seed pixel position
% CM: correlation matrix
% P: p-value

L=size(I,1);
C=size(I,2);
sp = squeeze(I(y,x,:));
I = reshape(I,L*C,size(I,3));
[CM,P] = corr(sp,I','Type','Pearson');
%[CM,P] = corr(sp,I','Type','Spearman');
CM = reshape(CM,L,C);
P = reshape(P,L,C);



% 


% VERSION 1
% function CM = OIA_SA_seedpixel(If,x,y)
% % function OIA_SA_seedpixel(I,x,y)
% % I: matrix of data
% % x,y: seed pixel position
% 
% ref = single(squeeze(If(y,x,:)));
% 
% CM = zeros(size(If,1),size(If,2));
% colormap jet
% for i = 1:size(If,1)
%     for j = 1:size(If,2)
%         s = single(squeeze(If(i,j,:)));
%         r = corrcoef(ref,s);
%         CM(i,j) = r(1,2);
%     end; 
%     imagesc(CM),
%     pause(.0001)
% end;
% %colorbar
% hold on, plot(x,y,'wo'), hold off
% title('seed pixel correlation')