%
% file to perform varimax rotation with Kaiser normalization
%
% function [AFIN]=varimkn(INPMAT)
function [AFIN T]=varimkn_beh(INPMAT)
% global pb;
% 
%   if ( exist('pb', 'var') )
%     pb.setStatus( 'Performing Kaiser Rotation' );
%     pb.refresh();
%   end;

  SUMI=sum((INPMAT.^2)');
  SUMINP=sqrt(SUMI);

  [m n]=size(INPMAT);
  for count1=1:m
    for count2=1:n
      KAISNOR(count1,count2)=INPMAT(count1,count2)./SUMINP(count1);
    end;
  end;

%   if ( exist('pb', 'var') )
%     pb.setStatus( 'Performing Varimax Rotation' );
%     pb.refresh();
%   end;

  [VARK T]=varim(KAISNOR);

  for count1=1:m
    for count2=1:n
      AFIN(count1,count2)=VARK(count1,count2).*SUMINP(count1);
    end;
  end;

