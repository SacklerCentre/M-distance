function [mdist_rS1, mdist_rS2] = get_mdist(X)

mdist_rS1 = (X.meta_ca - X.t2ca_rS1)/abs(X.meta_da);
mdist_rS2 = (-X.meta_ca + X.t2ca_rS2)/abs(X.meta_da);

end

