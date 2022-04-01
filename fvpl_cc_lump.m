function fit = fvpl_cc_lump(x)

% ps = length(x);
% fit = zeros(ps,1);
% for ii = 1:ps
%     f = analysisDriver('imex_param.in','responseIMEX.out',x(ii,:));
%     fit(ii) = f;
% end
outfile = [num2str(labindex) 'responseIMEX.out'];
fit = analysisDriver_cc_lump(outfile,x);
