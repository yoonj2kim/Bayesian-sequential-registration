function out = basis_bspline(f_domain, numBasis, df)
    out.x = f_domain;
    out.matrix = create_basismatrix(f_domain,numBasis,df);
end