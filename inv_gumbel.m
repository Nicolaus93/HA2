function gumbel=inv_gumbel(U,Beta,Mu)
%Inverse Gumbel Distribution used to similuate Gumbel samples

    gumbel=Mu-Beta*log(-log(U));