function varargout = shsemilogx(x,y,e,varargin)
%% H = shsemilogx(x,y,e,opt,varargin)
% This function creates a shaded error-bar plot on a log x axis. See help
% for `shplot` for more information.
%
% Written by Marcin Konowalczyk
% Timmel Group @ Oxford University

H = shplot(x,y,e,varargin{:});
set(gca,'XScale','log');
set(gca,'YScale','linear');
if nargout > 0; varargout = {H}; end
end