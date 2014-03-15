function y = up(x,M,FLAG)
% y = up(x,M); Upsample x by M
%
% Input:  x    - input signal
%         M    - upsampling factor
%         FLAG - 0/1, 0-dont append zeros, 1-append zeros so that
%                length(y)=length(x). Default FLAG=1
% Output: y - output signal
%
% Example:
%
% See also: down
%
%

if (nargin < 3)
  FLAG=1;
end;
y_len=length(x)*M;
if (FLAG==1)
y=zeros(1,y_len);
end;
y(1:M:y_len)=x;


