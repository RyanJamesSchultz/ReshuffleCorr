function [Tx,X,Xconf]=ReshuffleCorr(T, D1, D2, N, conf)
  % ReshuffleCorr(T, D1, D2, N, conf)
  %
  % This function performs the cross-correlation of two datasets with
  % reshuffling tests to determine confidence intervals.
  % Additional details are found within Schultz & Telesca [2017].
  %
  % Input data:
  % -T        -- time vector (shared between D1 & D2).
  % -D1       -- first data vector.
  % -D2       -- second data vector.
  % -count    -- number of surrogates computed during reshuffling test.
  % -conf     -- vector of confidence intervals to test (between 0-1).
  %
  % NOTES:
  %  - T is assumed to have uniform spacing.
  %  - D1 and D2 must share the same time axis, T.
  %  - This code assumes that D1 and D2 are vectors of real numbers.
  %  - negative lag means D1 arrives before D2.
  %  - confidence intervals require a sufficient number of surrogates.
  %
  % We kindly ask that if the user finds the use of this code beneficial 
  % that a citation to this code is included in their work:
  % Schultz, R., & Telesca, L. (2017). The Cross-Correlation and 
  % Reshuffling Tests in Discerning Induced Seismicity, [In preparation].
  %
  % See the above paper and references therein for examples of the
  % implementation of this code.
  %
  %
  % This program is free software: you can redistribute it and/or modify
  % it under the terms of the GNU General Public License as published
  % by the Free Software Foundation, either version 3 of the License, or
  % any later version.
  % This program is distributed in the hope that it will be useful,
  % but WITHOUT ANY WARRANTY; without even the implied warranty of
  % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  % GNU General Public License for more details:
  % http://www.gnu.org/licenses/
  
  
  % Check dimensions and sizes of input vectors.
  if(   (length(D1)~=length(D2)) || (length(D1)~=length(T)) )
      fprintf('ERROR:\n');
      fprintf('Improper input lengths\n');
      Tx=0; X=0; Xconf=0;
      return;
  end;
  
  % Check for 0-1 input on conf.
  if(any( (conf>1)|(conf<0)  ))
      fprintf('ERROR:\n');
      fprintf('Improper confidence intervals (must be between 0 and 1).\n');
      Tx=0; X=0; Xconf=0;
      return;
  end;
  
  % If they aren't already, change input vectors to column vectors.
  if(~iscolumn(T))
      T=T';
  end;
  if(~iscolumn(D1))
      D1=D1';
  end;
  if(~iscolumn(D2))
      D2=D2';
  end;
  
  % Calculate standardized time series and its cross-correlation.
  d1=D1-mean(D1);
  d2=D2-mean(D2);
  [X,Nx]=xcorr(d1,d2);
  X=X/( sqrt(xcorr(d1,d1,0))*sqrt(xcorr(d2,d2,0)) );
  
  % Determine output time axis.
  Num=max([length(D1),length(D2)]);
  dT=T(2)-T(1);
  Tx=Nx*dT;
  %num=floor(Num-1/2);
  %Tx=-num*dT:dT:dT*num;
  
  % Prep, calculate surogates, and then sort.
  XM=zeros(length(X), N);
  D_fft=abs(fft(d1));
  for i=1:N
      Di=d1(randperm(Num));
      Di=ifft(  D_fft.*exp(1i*angle(fft(Di))), 'symmetric');
      Xi=xcorr(Di,d2)/( sqrt(xcorr(Di,Di,0))*sqrt(xcorr(d2,d2,0)) );
      
      XM(:,i)=Xi;
  end;
  XM=sort(XM,2);
  
  % Determine confidence curves based on sorted surrogates and given confidence intervals.
  Xconf=zeros(length(X), length(conf));
  for i=1:length(conf)
      Xconf(:,i)=XM(:,floor(conf(i)*N));
  end;
  
return;
