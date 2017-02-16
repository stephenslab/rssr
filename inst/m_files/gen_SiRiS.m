function [SiRiS] = gen_SiRiS(R,se)
  Si 	  = 1 ./ se(:);
  [p,tp]=size(se);
  SiRiS   = full(sparse(repmat(Si, 1, p) .* R .* repmat(Si', p, 1))); 
end
