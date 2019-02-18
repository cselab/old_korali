function varargout = shsemilogy(x,y,e,varargin)
%% H = shsemilogy(x,y,e,opt,varargin)
% This function creates a shaded error-bar plot on a log y axis. See help
% for `shplot` for more information.
%
% Written by Marcin Konowalczyk
% Timmel Group @ Oxford University

H = shplot(x,y,e,varargin{:});
set(gca,'XScale','linear');
set(gca,'YScale','log');
if nargout > 0; varargout = {H}; end
end